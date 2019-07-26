//
// Created by Luis Yanes (EI) on 15/03/2018.
//

#pragma once

#include <cstdint>
#include <unordered_map>
#include <iostream>
#include <sdglib/factories/KMerIDXFactory.hpp>
#include <sdglib/utilities/omp_safe.hpp>
#include <sdglib/types/KmerTypes.hpp>
#include <sdglib/readers/FileReader.hpp>
#include <sdglib/factories/KmerPosFactory.hpp>
#include <sdglib/utilities/io_helpers.hpp>

class SequenceDistanceGraph;
class UniqueKmerIndex {
public:
    using Map = std::unordered_map<uint64_t, graphStrandPos>;
    using pair = std::pair<uint64_t, graphStrandPos>;
    using const_iterator = std::unordered_map<uint64_t, graphStrandPos>::const_iterator;
    UniqueKmerIndex(const SequenceDistanceGraph &sg, uint8_t k=31);

    const_iterator find(const uint64_t hash) const {
        return kmer_to_graphposition.find(hash);
    };

    const_iterator end() const {
        return kmer_to_graphposition.cend();
    };

    uint8_t get_k() const {
        return k;
    }

    std::tuple<bool, graphStrandPos> find_unique_kmer_in_graph(const uint64_t kmer) const {
        const auto nk = find(kmer);
        auto exists = end() != nk;
        graphStrandPos p;
        if (exists) {
            p = nk->second;
        }
        return std::make_tuple(exists, p);
    }

    bool is_unmappable(sgNodeID_t id) const {
        return 0 == unique_kmers_per_node[std::abs(id)];
    }

    uint64_t unique_kmers_in_node(const sgNodeID_t node) const {
        return unique_kmers_per_node[std::abs(node)];
    }

    uint64_t total_kmers_in_node(const sgNodeID_t node) const {
        return total_kmers_per_node[std::abs(node)];
    }

    std::string kmer_to_string(uint64_t kmer, int K) {
        static char nucleotides [4] = {'A', 'C', 'G', 'T'};
        std::string s;
        for (int i = 0; i < K; ++i) {
            s += nucleotides[kmer % 4];
            kmer = kmer / 4;
        }
        return s;
    }

    const Map& getMap() const {return kmer_to_graphposition; }
private:
    Map kmer_to_graphposition;
    uint8_t k = 31;
    std::vector<uint64_t> unique_kmers_per_node;
    std::vector<uint64_t> total_kmers_per_node;

};

class Unique63merIndex {
public:
    using Map = std::unordered_map<__uint128_t, graphStrandPos, int128_hash>;
    using pair = std::pair<__uint128_t, graphStrandPos>;
    using const_iterator = Map::const_iterator;

    Unique63merIndex(const SequenceDistanceGraph &sg, const uint8_t k=63);

    const_iterator find(const __uint128_t hash) const {
        return kmer_to_graphposition.find(hash);
    };

    const_iterator end() const {
        return kmer_to_graphposition.cend();
    };

    uint8_t get_k() const {
        return k;
    }

    std::tuple<bool, graphStrandPos> find_unique_kmer_in_graph(const __uint128_t kmer) const {
        const auto nk = find(kmer);
        auto exists = end() != nk;
        graphStrandPos p;
        if (exists) {
            p = nk->second;
        }
        return std::make_tuple(exists, p);
    }

    bool is_unmappable(sgNodeID_t id) const {
        return 0 == unique_kmers_per_node[std::abs(id)];
    }

    uint64_t unique_kmers_in_node(const sgNodeID_t node) const {
        return unique_kmers_per_node[std::abs(node)];
    }

    uint64_t total_kmers_in_node(const sgNodeID_t node) const {
        return total_kmers_per_node[std::abs(node)];
    }

    std::string kmer_to_string(uint64_t kmer, int K) {
        static char nucleotides [4] = {'A', 'C', 'G', 'T'};
        std::string s;
        for (int i = 0; i < K; ++i) {
            s += nucleotides[kmer % 4];
            kmer = kmer / 4;
        }
        return s;
    }

    //This should get deprecated ASAP
    const Map& getMap() const {return kmer_to_graphposition; }

private:
    Map kmer_to_graphposition;
    uint8_t k;
    std::vector<uint64_t> unique_kmers_per_node;
    std::vector<uint64_t> total_kmers_per_node;

};
