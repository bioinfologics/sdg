//
// Created by Bernardo Clavijo (EI) on 19/01/2018.
//

#ifndef BSG_GRAPHPARTITIONER_HPP
#define BSG_GRAPHPARTITIONER_HPP


#include "SequenceGraph.hpp"
#include "PairedReadMapper.hpp"
#include "KmerCompressionIndex.hpp"


class GraphPartitioner {
public:
    //instantiate receiving a graph and a 10x mapper, for now
    GraphPartitioner(SequenceGraph &_sg, std::vector<PairedReadMapper> & _rms, KmerCompressionIndex &_kci) : sg(_sg),rmappers(_rms),kci(_kci){};

    //gets all tag_patterns for a subgraph (does not return tag IDs, node order as per subgraph)
    std::vector<std::vector<bool>> tags_patterns(const SequenceSubGraph subgraph);

    //score_partition_set
    std::pair<float,float> score_partition_set(const SequenceSubGraph subgraph, const std::vector<SequenceSubGraph> partitions, const std::vector<std::vector<bool>> tag_patterns={});

    //generate partitions: receiving a list of nodes, create set list of partitions, return them as SequenceGraphPaths
    std::vector<SequenceSubGraph> generate_partitions(const SequenceSubGraph subgraph);

private:
    SequenceGraph &sg;
    std::vector<PairedReadMapper> &rmappers;
    KmerCompressionIndex &kci;


};


#endif //BSG_GRAPHPARTITIONER_HPP
