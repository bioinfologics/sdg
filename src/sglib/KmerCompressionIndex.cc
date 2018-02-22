//
// Created by Bernardo Clavijo (EI) on 03/12/2017.
//


#include "KmerCompressionIndex.hpp"
#include <atomic>


void KmerCompressionIndex::index_graph(){
    sglib::OutputLog(sglib::INFO) << "Indexing graph, Counting..."<<std::endl;
    const int k = 31;
    const int max_coverage = 1;
    uint64_t total_k=0;
    for (auto &n:sg.nodes) if (n.sequence.size()>=k) total_k+=n.sequence.size()+1-k;
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
        else if (*wi<*ri) {++wi; *wi=*ri;}
        else if (*wi==*ri){wi->merge(*ri);++ri;}
    }

    graph_kmers.resize(wi+1-graph_kmers.begin());
    sglib::OutputLog(sglib::INFO)<<graph_kmers.size()<<" kmers in index, "<< sizeof(KmerCount) << " bytes per kmer"<<std::endl;
    //TODO: remove kmers with more than X in count

//    std::vector<uint64_t> uniqKmer_statistics(kmerCount_SMR.summaryStatistics());
//    std::cout << "Number of " << int(k) << "-kmers seen in assembly " << uniqKmer_statistics[0] << std::endl;
//    std::cout << "Number of contigs from the assembly " << uniqKmer_statistics[2] << std::endl;
}

void KmerCompressionIndex::reindex_graph(){
    const int k = 31;
    const int max_coverage = 20;
    const std::string output_prefix("./");
    SMR<KmerCount,
    KmerCountFactory<FastaRecord>,
    GraphNodeReader<FastaRecord>,
    FastaRecord, GraphNodeReaderParams, KmerCountFactoryParams> kmerCount_SMR({1, sg}, {k}, max_mem, 0, max_coverage,
                                                                              output_prefix);



    std::cout << "Indexing graph kmers... " << std::endl;
    auto new_graph_kmers = kmerCount_SMR.process_from_memory();
    uint64_t deleted=0,changed=0,equal=0;
    //std::sort(new_graph_kmers.begin(),new_graph_kmers.end());
    for (auto i=0,j=0;i<graph_kmers.size() and j<new_graph_kmers.size();++j){
        while (i<graph_kmers.size() and graph_kmers[i].kmer<new_graph_kmers[j].kmer) {
            graph_kmers[i].count=0;
            ++deleted;
            ++i;
        }
        if (i<graph_kmers.size() and graph_kmers[i].kmer==new_graph_kmers[j].kmer){
            if (graph_kmers[i].count==new_graph_kmers[j].count) ++equal;
            else {
                graph_kmers[i].count=new_graph_kmers[j].count;
                ++changed;
            }
            ++i;
        }
    }
    std::cout << deleted << " deleted,   "<<changed<<" changed,   "<<equal<<" equal"<<std::endl;
    std::vector<uint64_t> uniqKmer_statistics(kmerCount_SMR.summaryStatistics());
    std::cout << "Number of " << int(k) << "-kmers seen in assembly " << uniqKmer_statistics[0] << std::endl;
    std::cout << "Number of contigs from the assembly " << uniqKmer_statistics[2] << std::endl;
}

void KmerCompressionIndex::load_from_disk(std::string filename) {
    std::ifstream inf(filename);
    //read-to-tag
    uint64_t kcount;
    inf.read(( char *) &kcount,sizeof(kcount));
    graph_kmers.resize(kcount);
    inf.read(( char *) graph_kmers.data(),sizeof(KmerCount)*kcount);
    //read-to-node
    uint64_t ccount;
    inf.read(( char *) &ccount,sizeof(ccount));
    for (auto i=0;i<ccount;++i) {
        read_counts.emplace_back();
        read_counts.back().resize(kcount);
        inf.read(( char *) read_counts.back().data(), sizeof(uint16_t) * kcount);
    }

}

void KmerCompressionIndex::save_to_disk(std::string filename) {
    std::ofstream of(filename);
    //read-to-tag
    uint64_t kcount=graph_kmers.size();
    of.write((const char *) &kcount,sizeof(kcount));
    of.write((const char *) graph_kmers.data(),sizeof(KmerCount)*kcount);
    //read-to-node
    uint64_t ccount=read_counts.size();
    of.write((const char *) &ccount,sizeof(ccount));
    for (auto i=0;i<ccount;++i) {
        of.write((const char *) read_counts[i].data(), sizeof(uint16_t) * kcount);
    }
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
    std::cout << "Mean coverage for unique kmers:   " << ((double)tuniq)/cuniq <<std::endl;
    std::cout << "Median coverage for unique kmers: " << median <<std::endl;
    std::cout << "Mode coverage for unique kmers:   " << mode <<std::endl;

    if (median<.9*mode or median>.9*mode ) std::cout<<"WARNING -> median and mode highly divergent"<<std::endl;
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

double KmerCompressionIndex::compute_compression_for_node(sgNodeID_t _node, uint16_t max_graph_freq) {

    auto & node=sg.nodes[_node>0 ? _node:-_node];

    std::vector<uint64_t> nkmers;
    StringKMerFactory skf(node.sequence,31);
    skf.create_kmers(nkmers);

    uint64_t kcount=0,kcov=0;
    for (auto &kmer : nkmers){
        auto nk = std::lower_bound(graph_kmers.begin(), graph_kmers.end(), KmerCount(kmer,0));
        if (nk!=graph_kmers.end() and nk->kmer == kmer and nk->count==1) {
            ++kcount;
            kcov+=read_counts[0][nk-graph_kmers.begin()];
        }
    }

    return (((double) kcov)/kcount )/uniq_mode;
}
