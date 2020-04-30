//
// Created by Bernardo Clavijo (EI) on 30/04/2020.
//

#pragma once
#include <sdglib/graph/DistanceGraph.hpp>

class PerfectMatchPart{
public:
    void extend(const std::string & readseq,const std::string & nodeseq);


    sgNodeID_t node;

    uint64_t offset;//will only be set if first part
    int64_t previous_part;//will only be set if not first part
    uint64_t read_position;//position of the last matched base in read
    uint64_t node_position;//postiion of the last matched base in node (canonical orientation)

    bool completed_node=false;
    bool completed_read=false;
    bool extended=false;
    bool invalid=false;
};

class PerfectMatchExtender{
public:
    PerfectMatchExtender(DistanceGraph & _dg, uint8_t _k):dg(_dg),k(_k){};

    void reset(std::string _readseq);
    void add_starting_match(sgNodeID_t node_id, uint64_t read_offset, uint64_t node_offset);
    void extend_fw();
    std::vector<sgNodeID_t> best_path(); //Todo: return pointer to the last part?



    DistanceGraph & dg;
    uint8_t k;
    std::vector<PerfectMatchPart> matchparts;
    uint64_t last_readpos;
    uint64_t last_nodepos;
    std::string readseq;

};

class PerfectMatcher {
    //Create with a graph and parameters for the index.

    //set a sequence to map

    //get_next_match returns a PerfectMatchPath object that points to the pme and enables to get the match parts

};