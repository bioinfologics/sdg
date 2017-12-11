//
// Created by Bernardo Clavijo (EI) on 03/12/2017.
//


#include "KmerCompressionIndex.hpp"



KmerCompressionIndex::KmerCompressionIndex(SequenceGraph &_sg, uint64_t _max_mem): sg(_sg) {
    max_mem = _max_mem;
}

void KmerCompressionIndex::index_graph(){
    const int k = 31;
    const int max_coverage = 1;
    const std::string output_prefix("./");
    SMR<KmerCount,
    KmerCountFactory<FastaRecord>,
    GraphNodeReader<FastaRecord>,
    FastaRecord, GraphNodeReaderParams, KmerCountFactoryParams> kmerCount_SMR({1, sg}, {k}, max_mem, 0, max_coverage,
                                                                          output_prefix);



    std::cout << "Indexing assembly... " << std::endl;
    graph_kmers = kmerCount_SMR.process_from_memory();

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

void KmerCompressionIndex::add_counts_from_file(std::string filename) {

    std::cout<<"counting from file"<<std::endl;

    FastqReader<FastqRecord> fastqReader({0},filename);
    FastqRecord read;
    std::vector<KmerCount> readkmers;
    KmerCountFactory<FastqRecord> kf({31});
    uint64_t present=0,absent=0,rp=0;
    bool c ;
    c = fastqReader.next_record(read);
    while (c) {
        readkmers.clear();
        //process tag if 10x! this way even ummaped reads get tags
        kf.setFileRecord(read);
        kf.next_element(readkmers);

        for (auto &rk:readkmers) {
            auto nk = std::lower_bound(graph_kmers.begin(), graph_kmers.end(), rk);
            if (nk!=graph_kmers.end() and nk->kmer == rk.kmer) {
                auto offset = nk - graph_kmers.begin();
                read_counts.back()[offset]++;
                ++present;
            }
            else ++absent;
        }
        ++rp;
        if (rp % 100000 == 0) std::cout << rp << " reads processed "<< present <<" / " << present+absent << " kmers found" << std::endl;
        c = fastqReader.next_record(read);
    }
    // somehow for my test data with 700 reads, totak count is 834 for r2...
    std::cout << rp << " reads processed "<< present <<" / " << present+absent << " kmers found" << std::endl;
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