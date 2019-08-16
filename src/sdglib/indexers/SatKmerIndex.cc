//
// Created by Luis Yanes (EI) on 2019-06-25.
//

#include "SatKmerIndex.hpp"
#include <sdglib/graph/SequenceDistanceGraph.hpp>

SatKmerIndex::SatKmerIndex(const SequenceDistanceGraph &sg, uint8_t k, uint8_t filter_limit)  : k(k) {
    if (k > 15) {
        throw std::runtime_error(
                "You are trying to use K>15, which is not supported by this structure. "
                "Please consider NKmerIndex or UniqueKmerIndex as alternatives");
    }
    // this can be parallelised by contig by aggregating different k_usage vectors on the first step, second and third steps are trickier
    // This two-pass approach could be easily adapted to a single memory block rather than vector of vectors
    uint64_t filtered_kmers(0),total_kmers(0);
    //---- First Step, count k-mer occupancy ----//
    kmerEnd.resize(std::pow(4,k));
    std::vector<kmerPos> all_kmers;
    all_kmers.reserve(100000000);
    std::vector<std::pair<bool,uint64_t > > contig_kmers;
    contig_kmers.reserve(1000000); //1Mbp per contig to start with?
    StringKMerFactory skf(k);
    for (sgNodeID_t n = 1; n < sg.nodes.size(); ++n) {
        if (sg.nodes[n].sequence.size() >= k) {
            contig_kmers.clear();
            skf.create_kmers(sg.nodes[n].sequence, contig_kmers);
            int k_i(0);
            for (const auto &kmer:contig_kmers) {
                all_kmers.emplace_back(kmer.second, n, kmer.first ? k_i + 1 : -(k_i + 1));
                k_i++;
            }
            total_kmers+=contig_kmers.size();
        }
    }
    //---- Second Step, reserve space for each vector in structure (avoiding reallocations and such)----//
    sdglib::sort(all_kmers.begin(), all_kmers.end(), kmerPos::byKmerContigOffset());
    auto witr = all_kmers.begin();
    auto ritr = witr;
    uint64_t ckmer(0);
    for (; ritr != all_kmers.end();) {

        while (ckmer != ritr->kmer) {
            ++ckmer;
            kmerEnd[ckmer]=kmerEnd[ckmer-1];
        }
        auto bitr = ritr;
        while (ritr != all_kmers.end() and bitr->kmer == ritr->kmer) {
            ++ritr;
        }
        if (ritr - bitr < filter_limit) {
            filtered_kmers+=1;
            kmerEnd[ckmer] += ritr-bitr;
            while (bitr != ritr) {
                contig_offsets.emplace_back(bitr->contigID, bitr->offset, sg.nodes[llabs(bitr->contigID)].sequence.size() - std::abs(bitr->offset));
                ++witr;
                ++bitr;
            }
        }
    }
    while (ckmer < std::pow(4,k)) {
        ++ckmer;
        kmerEnd[ckmer] = kmerEnd[ckmer-1];
    }
    assert(kmerEnd[std::pow(4,k)-1] == contig_offsets.size());
    std::vector<kmerPos>().swap(all_kmers);
}