//
// Created by Bernardo Clavijo (EI) on 10/05/2018.
//

#ifndef BSG_READMAPPING_HPP
#define BSG_READMAPPING_HPP
#include "sglib/SequenceGraph.hpp"



class ReadMapping {
public:
    ReadMapping(){
        //just clean the structure, OSX doesn't give you clean memory
        bzero(this, sizeof(ReadMapping));
    }
    bool operator==(const ReadMapping &other){
        return this==&other;
    };
    bool operator<(const ReadMapping &other) const {
        if (node!=other.node) return node<other.node;
        return read_id<other.read_id;
    };
    void merge(const ReadMapping &other){};

    sgNodeID_t node;
    uint64_t read_id;
    int32_t first_pos;
    int32_t last_pos;
    int32_t unique_matches;
    bool rev=false;

};
#endif //BSG_READMAPPING_HPP
