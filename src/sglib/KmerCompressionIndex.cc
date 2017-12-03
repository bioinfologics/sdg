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

    std::vector<KmerCount> unique_kmers;

    std::cout << "Indexing assembly... " << std::endl;
    unique_kmers = kmerCount_SMR.read_from_file("kcompidx");

    std::vector<uint64_t> uniqKmer_statistics(kmerCount_SMR.summaryStatistics());
    std::cout << "Number of " << int(k) << "-kmers seen in assembly " << uniqKmer_statistics[0] << std::endl;
    std::cout << "Number of contigs from the assembly " << uniqKmer_statistics[2] << std::endl;
}

void KmerCompressionIndex::add_counts_from_file(std::string filename) {

}