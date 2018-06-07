//
// Created by Ben Ward (EI) on 07/02/2018.
//

#ifndef BSG_UNIQUEKMERINDEX_H
#define BSG_UNIQUEKMERINDEX_H

#include <sglib/types/KmerTypes.hpp>
#include "sglib/graph/SequenceGraph.hpp"

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

    bool traverse_dark_nodes(const sgNodeID_t seed, const sgNodeID_t goal);
};



#endif //BSG_UNIQUEKMERINDEX_H
