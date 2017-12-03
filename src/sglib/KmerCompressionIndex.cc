//
// Created by Bernardo Clavijo (EI) on 03/12/2017.
//


#include "KmerCompressionIndex.hpp"



KmerCompressionIndex::KmerCompressionIndex(SequenceGraph &_sg, uint64_t _max_mem): sg(_sg) {
    max_mem=_max_mem;

    const int k = 31;
    const int max_coverage = 1;
    const std::string output_prefix("./");
    SMR<KmerCount,
    KmerCountFactory<FastaRecord>,
    GraphNodeReader<FastaRecord>,
    FastaRecord, GraphNodeReaderParams, KmerCountFactoryParams> kmerCount_SMR({1, sg}, {k}, max_mem, 0, max_coverage,
                                                                          output_prefix);



    std::cout << "Indexing assembly... " << std::endl;
    graph_kmers = kmerCount_SMR.read_from_file("kcompidx");

    std::vector<uint64_t> uniqKmer_statistics(kmerCount_SMR.summaryStatistics());
    std::cout << "Number of " << int(k) << "-kmers seen in assembly " << uniqKmer_statistics[0] << std::endl;
    std::cout << "Number of contigs from the assembly " << uniqKmer_statistics[2] << std::endl;
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