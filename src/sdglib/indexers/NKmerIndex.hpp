//
// Created by Luis Yanes (EI) on 17/08/2018.
//

#pragma once

#include <vector>
#include <sdglib/utilities/omp_safe.hpp>
#include <sdglib/types/KmerTypes.hpp>
#include <sdglib/factories/KMerFactory.hpp>
#include <sdglib/utilities/OutputLog.hpp>
#include <sdglib/graph/SequenceDistanceGraph.hpp>
#include <sdglib/bloom/BloomFilter.hpp>

class SequenceDistanceGraph;

class NKmerIndex {
    BloomFilter bfilter;
    std::vector<kmerPos> assembly_kmers;
    uint8_t k=0;
    const SequenceDistanceGraph &sg;
public:
    using const_iterator = std::vector<kmerPos>::const_iterator;

    explicit NKmerIndex(const SequenceDistanceGraph &_sg,uint8_t k=15, int filter_limit = 200);

    bool empty() const { return assembly_kmers.empty(); }
    const_iterator begin() const {return assembly_kmers.cbegin();}
    const_iterator end() const {return assembly_kmers.cend();}

    bool filter(const uint64_t kmer) const {
        return bfilter.contains(kmer);
    }

    const_iterator find(const uint64_t kmer) const {
        if (bfilter.contains(kmer)) { return std::lower_bound(assembly_kmers.cbegin(), assembly_kmers.cend(), kmer); }
        return assembly_kmers.cend();
    }

    void get_all_kmer_matches(std::vector<std::vector<std::pair<int32_t, int32_t>>> & matches, std::vector<std::pair<bool, uint64_t>> & seq_kmers);

};


class NKmerIndex128 {
    BloomFilter bfilter;
    std::vector<kmerPos128> assembly_kmers;
    uint8_t k=0;
    const SequenceDistanceGraph &sg;
public:
    using const_iterator = std::vector<kmerPos128>::const_iterator;

    explicit NKmerIndex128(const SequenceDistanceGraph &_sg,uint8_t k=63, int filter_limit = 200);

    bool empty() const { return assembly_kmers.empty(); }
    const_iterator begin() const {return assembly_kmers.cbegin();}
    const_iterator end() const {return assembly_kmers.cend();}

    bool filter(const __uint128_t kmer) const {
        return bfilter.contains(kmer);
    }

    const_iterator find(const __uint128_t kmer) const {
        if (bfilter.contains(kmer)) { return std::lower_bound(assembly_kmers.cbegin(), assembly_kmers.cend(), kmer); }
        return assembly_kmers.cend();
    }

    void get_all_kmer_matches(std::vector<std::vector<std::pair<int32_t, int32_t>>> & matches, std::vector<std::pair<bool, __uint128_t>> & seq_kmers);

};
