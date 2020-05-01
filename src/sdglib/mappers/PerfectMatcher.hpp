//
// Created by Bernardo Clavijo (EI) on 30/04/2020.
//

#pragma once
#include <sdglib/graph/DistanceGraph.hpp>
#include <sdglib/indexers/NKmerIndex.hpp>
#include <memory>

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
    void set_best_path(); //Todo: return pointer to the last part?



    DistanceGraph & dg;
    uint8_t k;
    std::vector<PerfectMatchPart> matchparts;
    std::vector<sgNodeID_t> best_path;
    uint64_t last_readpos;
    uint64_t last_nodepos;
    std::string readseq;

};

/*class PerfectMatcher {
public:
    //Create with a graph and parameters for the index.
    PerfectMatcher(DistanceGraph &_dg,std::shared_ptr<NKmerIndex> _nki):dg(_dg),nki(_nki){};

    //Create with a graph and a pointer to the index.
    PerfectMatcher(DistanceGraph &_dg,uint8_t _k, uint16_t _max_freq):dg(_dg),nki(std::make_shared<NKmerIndex>(dg.sdg,_k,_max_freq)){};

    //set a sequence to map
    void set_sequence();

    //get_next_match returns a PerfectMatchPath object that points to the pme and enables to get the match parts
    PerfectMatchPart & get_next_part();

    void reset_part_getter()

    std::shared_ptr<NKmerIndex> nki;
    DistanceGraph & dg;
};*/