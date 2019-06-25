//
// Created by Luis Yanes (EI) on 17/08/2018.
//

#ifndef BSG_NKMERINDEX_HPP
#define BSG_NKMERINDEX_HPP

#include <vector>
#include <sdglib/utilities/omp_safe.hpp>
#include <sdglib/types/KmerTypes.hpp>
#include <sdglib/factories/KMerFactory.hpp>
#include <sdglib/logger/OutputLog.hpp>
#include <sdglib/graph/SequenceDistanceGraph.hpp>
#include <sdglib/bloom/BloomFilter.hpp>

class NKmerIndex {
    BloomFilter filter;
    std::vector<kmerPos> assembly_kmers;
    uint8_t k=15;
public:
    using const_iterator = std::vector<kmerPos>::const_iterator;

    NKmerIndex(){}
    explicit NKmerIndex(uint8_t k) : k(k), filter(70*1024*1024) {}

    uint64_t filter_kmers(std::vector<kmerPos> &kmers, int max_kmer_repeat) {
        uint64_t total_kmers(0);
        auto witr = kmers.begin();
        auto ritr = witr;
        for (; ritr != kmers.end();) {
            auto bitr = ritr;
            while (ritr != kmers.end() and bitr->kmer == ritr->kmer) {
                ++ritr;
            }
            if (ritr - bitr < max_kmer_repeat) {
                total_kmers+=1;
                while (bitr != ritr) {
                    *witr = *bitr;
                    ++witr;
                    ++bitr;
                }
            }
        }
        kmers.resize(witr-kmers.begin());
        return total_kmers;
    }

    void generate_index(const SequenceDistanceGraph &sg, int filter_limit = 200, bool verbose=true);

    bool empty() const { return assembly_kmers.empty(); }
    const_iterator begin() const {return assembly_kmers.cbegin();}
    const_iterator end() const {return assembly_kmers.cend();}

    const_iterator find(const uint64_t kmer) const {
        if (filter.contains(kmer)) { return std::lower_bound(assembly_kmers.cbegin(), assembly_kmers.cend(), kmer); }
        return assembly_kmers.cend();
//        return std::lower_bound(assembly_kmers.cbegin(), assembly_kmers.cend(), kmer);
    }
};


#endif //BSG_NKMERINDEX_HPP
