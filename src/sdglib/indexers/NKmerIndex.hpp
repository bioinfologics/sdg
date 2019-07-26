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
    BloomFilter filter;
    std::vector<kmerPos> assembly_kmers;
    uint8_t k=0;
public:
    using const_iterator = std::vector<kmerPos>::const_iterator;

    NKmerIndex(){};
    NKmerIndex(const SequenceDistanceGraph &sg,uint8_t k=15, int filter_limit = 200);

    bool empty() const { return assembly_kmers.empty(); }
    const_iterator begin() const {return assembly_kmers.cbegin();}
    const_iterator end() const {return assembly_kmers.cend();}

    const_iterator find(const uint64_t kmer) const {
        if (filter.contains(kmer)) { return std::lower_bound(assembly_kmers.cbegin(), assembly_kmers.cend(), kmer); }
        return assembly_kmers.cend();
    }
};
