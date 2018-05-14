//
// Created by Bernardo Clavijo (EI) on 10/05/2018.
//

#ifndef BSG_READMAPPER_HPP
#define BSG_READMAPPER_HPP
#include "sglib/SequenceGraph.hpp"



class ReadMapper {
public:
    ReadMapper(){
        //just clean the structure, OSX doesn't give you clean memory
        bzero(this, sizeof(ReadMapper));
    }
    bool operator==(const ReadMapper &other){
        return this==&other;
    };
    bool operator<(const ReadMapper &other) const {
        if (node!=other.node) return node<other.node;
        return read_id<other.read_id;
    };
    void merge(const ReadMapper &other){};

    sgNodeID_t node;
    uint64_t read_id;
    int32_t first_pos;
    int32_t last_pos;
    int32_t unique_matches;
    bool rev=false;

};
#endif //BSG_READMAPPER_HPP
