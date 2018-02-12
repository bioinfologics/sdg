//
// Created by Ben Ward (EI) on 07/02/2018.
//

#ifndef BSG_UNIQUEKMERINDEX_H
#define BSG_UNIQUEKMERINDEX_H

#include "SequenceGraph.hpp"

struct graphPosition{
    sgNodeID_t node;
    uint32_t pos;
};

class UniqueKmerIndex {
private:
    const SequenceGraph& sg;
    std::unordered_map<uint64_t, graphPosition> kmer_to_node_map;
    std::vector<uint64_t> unique_kmers_per_node;
    std::vector<uint64_t> total_kmers_per_node;

public:
    explicit UniqueKmerIndex(const SequenceGraph& sg, uint8_t k);

    uint64_t unique_kmers_in_node(const sgNodeID_t node) const {
        return unique_kmers_per_node[std::abs(node)];
    }

    uint64_t total_kmers_in_node(const sgNodeID_t node) const {
        return total_kmers_per_node[std::abs(node)];
    }

    std::tuple<bool, graphPosition> find_unique_kmer_in_graph(uint64_t kmer) const;
};



#endif //BSG_UNIQUEKMERINDEX_H
