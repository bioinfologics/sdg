//
// Created by Bernardo Clavijo (EI) on 03/12/2017.
//


#include "KmerCompressionIndex.hpp"
#include <atomic>
#include <sglib/readers/FileReader.h>


void KmerCompressionIndex::index_graph(){
    sglib::OutputLog(sglib::INFO) << "Indexing graph, Counting..."<<std::endl;
    const int k = 31;
    uint64_t total_k=0;
    for (auto &n:sg.nodes) if (n.sequence.size()>=k) total_k+=n.sequence.size()+1-k;
    graph_kmers.clear();
    graph_kmers.reserve(total_k);
    FastaRecord r;
    KmerCountFactory<FastaRecord>kcf({k});
    for (sgNodeID_t n=1;n<sg.nodes.size();++n){
        if (sg.nodes[n].sequence.size()>=k){
            r.id=n;
            r.seq=sg.nodes[n].sequence;
            kcf.setFileRecord(r);
            kcf.next_element(graph_kmers);
        }
    }
    sglib::OutputLog(sglib::INFO)<<graph_kmers.size()<<" kmers in total"<<std::endl;
    sglib::OutputLog(sglib::INFO) << "  Sorting..."<<std::endl;
    std::sort(graph_kmers.begin(),graph_kmers.end());
    sglib::OutputLog(sglib::INFO) << "  Merging..."<<std::endl;
    auto wi=graph_kmers.begin();
    auto ri=graph_kmers.begin();
    while (ri<graph_kmers.end()){
        if (wi.base()==ri.base()) ++ri;
        else if (*wi<*ri) {++wi; *wi=*ri;++ri;}
        else if (*wi==*ri){wi->merge(*ri);++ri;}
    }

    graph_kmers.resize(wi+1-graph_kmers.begin());
    sglib::OutputLog(sglib::INFO)<<graph_kmers.size()<<" kmers in index"<<std::endl;
    //TODO: remove kmers with more than X in count

//    std::vector<uint64_t> uniqKmer_statistics(kmerCount_SMR.summaryStatistics());
//    std::cout << "Number of " << int(k) << "-kmers seen in assembly " << uniqKmer_statistics[0] << std::endl;
//    std::cout << "Number of contigs from the assembly " << uniqKmer_statistics[2] << std::endl;
}

void KmerCompressionIndex::reindex_graph(){
    sglib::OutputLog(sglib::INFO) << "Re-indexing graph, Counting..."<<std::endl;
    std::vector<uint8_t> new_counts(graph_kmers.size());

    const int k = 31;
    std::vector<KmerCount> nodekmers;
    FastaRecord r;
    KmerCountFactory<FastaRecord>kcf({k});
    for (sgNodeID_t n=1;n<sg.nodes.size();++n){
        if (sg.nodes[n].sequence.size()>=k){
            r.id=n;
            r.seq=sg.nodes[n].sequence;
            kcf.setFileRecord(r);
            kcf.next_element(nodekmers);
            for (auto &nk:nodekmers){
                auto gk=std::lower_bound(graph_kmers.begin(),graph_kmers.end(),nk);
                if (gk->kmer==nk.kmer and new_counts[gk-graph_kmers.begin()]<255) ++new_counts[gk-graph_kmers.begin()];
            }
            nodekmers.clear();
        }
    }
    sglib::OutputLog(sglib::INFO) << "  Updating counts..."<<std::endl;
    for (auto i=0;i<graph_kmers.size();++i){
        graph_kmers[i].count=new_counts[i];
    }
    sglib::OutputLog(sglib::INFO) << "Re-indexing done."<<std::endl;
}

void KmerCompressionIndex::read(std::ifstream &input_file) {
    uint64_t kcount;
    input_file.read(( char *) &kcount,sizeof(kcount));
    graph_kmers.resize(kcount);
    input_file.read(( char *) graph_kmers.data(),sizeof(KmerCount)*kcount);
    //read-to-node
    uint64_t ccount;
    input_file.read(( char *) &ccount,sizeof(ccount));
    for (auto i=0;i<ccount;++i) {
        read_counts.emplace_back();
        read_counts.back().resize(kcount);
        input_file.read(( char *) read_counts.back().data(), sizeof(uint16_t) * kcount);
    }
}

void KmerCompressionIndex::load_from_disk(std::string filename) {
    std::ifstream inf(filename);
    //read-to-tag
    read(inf);

}
void KmerCompressionIndex::write(std::ofstream &output_file) {
    uint64_t kcount=graph_kmers.size();
    output_file.write((const char *) &kcount,sizeof(kcount));
    output_file.write((const char *) graph_kmers.data(),sizeof(KmerCount)*kcount);
    uint64_t ccount=read_counts.size();
    output_file.write((const char *) &ccount,sizeof(ccount));
    for (auto i=0;i<ccount;++i) {
        output_file.write((const char *) read_counts[i].data(), sizeof(uint16_t) * kcount);
    }
}
void KmerCompressionIndex::save_to_disk(std::string filename) {
    std::ofstream of(filename);
    write(of);
}

void KmerCompressionIndex::start_new_count(){
    read_counts.emplace_back();
    read_counts.back().resize(graph_kmers.size(),0);
}

void KmerCompressionIndex::add_counts_from_file(std::vector<std::string> filenames) {


    uint64_t present(0), absent(0), rp(0);
    sglib::OutputLog(sglib::INFO)<<"Populating lookup map"<<std::endl;
    std::unordered_map<uint64_t,uint64_t> kmer_map;
    kmer_map.reserve(graph_kmers.size());
    for (uint64_t i=0;i<graph_kmers.size();++i) kmer_map[graph_kmers[i].kmer]=i;
    sglib::OutputLog(sglib::INFO)<<"Map populated, processing files"<<std::endl;
    for (auto filename:filenames) {
        sglib::OutputLog(sglib::INFO) << "Counting from file: " << filename << std::endl;
        FastqReader<FastqRecord> fastqReader({0}, filename);

#pragma omp parallel shared(fastqReader)
        {
            uint64_t thread_present(0), thread_absent(0), thread_rp(0);
            const size_t local_kmers_size = 2000000;
            std::vector<uint64_t> found_kmers;
            found_kmers.reserve(local_kmers_size);
            FastqRecord read;
            std::vector<KmerCount> readkmers;
            KmerCountFactory<FastqRecord> kf({31});

            bool c;
#pragma omp critical(fastqreader)
            c = fastqReader.next_record(read);
            while (c) {
                readkmers.clear();
                kf.setFileRecord(read);
                kf.next_element(readkmers);

                for (auto &rk:readkmers) {
                    auto findk = kmer_map.find(rk.kmer);
                    if (kmer_map.end() != findk) {
                        //++thread_counts[findk->second];
                        found_kmers.emplace_back(findk->second);
                        if (found_kmers.size() == local_kmers_size) {
                            bool printstatus = false;
#pragma omp critical(results_merge)
                            {
                                auto &arc = read_counts.back();
                                for (auto &x:found_kmers) if (arc[x] < UINT16_MAX) ++arc[x];
                                if (rp / 1000000 != (rp + thread_rp) / 1000000) printstatus = true;
                                rp += thread_rp;
                                present += thread_present;
                                absent += thread_absent;
                            }
                            found_kmers.clear();
                            thread_absent = 0;
                            thread_present = 0;
                            thread_rp = 0;
                            if (printstatus)
                                sglib::OutputLog(sglib::INFO) << rp << " reads processed " << present << " / "
                                                              << present + absent << " kmers found" << std::endl;
                        }
                        ++thread_present;
                    } else ++thread_absent;


                }
                ++thread_rp;
#pragma omp critical(fastqreader)
                c = fastqReader.next_record(read);
            }
        }
        sglib::OutputLog(sglib::INFO) << rp << " reads processed " << present << " / " << present + absent
                                      << " kmers found" << std::endl;
    }
}



void KmerCompressionIndex::compute_compression_stats() {
    //compute mean, median and mode, as of now, only use the first read count
    //uint64_t graphcov[10]={0,0,0,0,0,0,0,0,0,0};
    //std::cout << "Coverage in graph:" <<std::endl;
    //for (auto &gk:graph_kmers) ++graphcov[(gk.count<10?gk.count-1:9)];
    //for (auto i=1;i<10;++i) std::cout << i <<":   "<<graphcov[i-1]<<std::endl;
    //std::cout <<"10+: "<<graphcov[9]<<std::endl;
    uint64_t covuniq[1001];
    for (auto &c:covuniq)c=0;
    uint64_t tuniq=0,cuniq=0;
    for (uint64_t i=0; i<graph_kmers.size(); ++i){
        if (graph_kmers[i].count==1){
            tuniq+=read_counts[0][i];
            ++cuniq;
            ++covuniq[(read_counts[0][i]<1000 ? read_counts[0][i] : 1000 )];
        }
    }
    uint64_t cseen=0,median=0;
    while (cseen<cuniq/2) {cseen+=covuniq[median];++median;};
    uint64_t mode=0;
    for (auto i=0;i<1000;++i) if (covuniq[i]>covuniq[mode]) mode=i;
    sglib::OutputLog()<<"KCI Mean coverage for unique kmers:   " << ((double)tuniq)/cuniq <<std::endl;
    sglib::OutputLog()<<"KCI Median coverage for unique kmers: " << median <<std::endl;
    sglib::OutputLog()<<"KCI Mode coverage for unique kmers:   " << mode <<std::endl;

    if (median<.9*mode or median>.9*mode ) sglib::OutputLog()<<"WARNING -> median and mode highly divergent"<<std::endl;
    uniq_mode=mode;

}

void KmerCompressionIndex::dump_histogram(std::string filename) {
    std::ofstream kchf(filename);
    uint64_t covuniq[1001];
    for (auto &c:covuniq)c=0;
    uint64_t tuniq=0,cuniq=0;
    for (uint64_t i=0; i<graph_kmers.size(); ++i){
        if (graph_kmers[i].count==1){
            tuniq+=read_counts[0][i];
            ++cuniq;
            ++covuniq[(read_counts[0][i]<1000 ? read_counts[0][i] : 1000 )];
        }
    }
    for (auto i=0;i<1000;++i) kchf<<i<<","<<covuniq[i]<<std::endl;
}

double KmerCompressionIndex::compute_compression_for_node(sgNodeID_t _node, uint16_t max_graph_freq, uint16_t dataset) {
    const int k=31;
    auto n=_node>0 ? _node:-_node;
    auto & node=sg.nodes[n];

    //eliminate "overlapping" kmers
    int32_t max_bw_ovlp=0;
    int32_t max_fw_ovlp=0;
    for (auto bl:sg.get_bw_links(n)) {
        if (bl.dist<0){
            auto ovl=-bl.dist+1-k;
            if (ovl>max_bw_ovlp) max_bw_ovlp=ovl;
        }
    }
    for (auto fl:sg.get_fw_links(n)) {
        if (fl.dist<0){
            auto ovl=-fl.dist+1-k;
            if (ovl>max_fw_ovlp) max_fw_ovlp=ovl;
        }
    }
    int64_t newsize=node.sequence.size();
    newsize=newsize-max_bw_ovlp-max_fw_ovlp;
    //if (n/10==50400){
    //    std::cout<<"node "<<n<<" size="<<node.sequence.size()<<" max_bw_olv="<<max_bw_ovlp<<" max_fw_ovl="<<max_fw_ovlp<<" newlength="<<newsize<<std::endl;
    //}
    if (newsize<k) return ((double)0/0);
    auto s=node.sequence.substr(max_bw_ovlp,newsize);
    std::vector<uint64_t> nkmers;
    StringKMerFactory skf(s,k);


    skf.create_kmers(nkmers);

    uint64_t kcount=0,kcov=0;
    for (auto &kmer : nkmers){
        auto nk = std::lower_bound(graph_kmers.begin(), graph_kmers.end(), KmerCount(kmer,0));
        if (nk!=graph_kmers.end() and nk->kmer == kmer and nk->count<=max_graph_freq) {
            kcount+=nk->count;
            kcov+=read_counts[dataset][nk-graph_kmers.begin()];
        }
    }

    return (((double) kcov)/kcount )/uniq_mode;
}

void KmerCompressionIndex::compute_all_nodes_kci(uint16_t max_graph_freq) {
    nodes_depth.resize(sg.nodes.size());
    for (auto n=1;n<sg.nodes.size();++n) {
        nodes_depth[n]=compute_compression_for_node(n, max_graph_freq);
    }
}