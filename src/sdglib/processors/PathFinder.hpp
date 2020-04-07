//
// Created by Bernardo Clavijo (EI) on 07/04/2020.
//

#pragma once


#include <sdglib/workspace/WorkSpace.hpp>

/*
 * This tries to choose the right path between two haplotype specific nodes, by using whatever data we can put to it.
 */

class PathFinder {
public:
    PathFinder(WorkSpace &_ws, sgNodeID_t _n1, sgNodeID_t _n2, std::vector<SequenceDistanceGraphPath> _paths, uint8_t _k):
        ws(_ws),n1(_n1),n2(_n2),paths(_paths),k(_k){
        index_paths();
    };
    void index_paths();
    std::vector<std::vector<uint64_t>> seq_to_pathpos(uint16_t path_id, std::string seq);
    std::unordered_map<uint64_t,std::vector<std::pair<uint16_t ,uint64_t >>> kmerpos; //kmer -> [(path_id, pos)]
    WorkSpace & ws;
    sgNodeID_t n1,n2;
    std::vector<SequenceDistanceGraphPath> paths;
    uint8_t k;
};


