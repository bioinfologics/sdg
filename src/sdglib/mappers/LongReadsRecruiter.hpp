//
// Created by Bernardo Clavijo (EI) on 16/01/2020.
//

#ifndef SDG_LONGREADSRECRUITER_HPP
#define SDG_LONGREADSRECRUITER_HPP


#include <sdglib/graph/SequenceDistanceGraph.hpp>
#include <sdglib/datastores/LongReadsDatastore.hpp>

//TODO: OMP
//TODO: dump/load
//TODO: translate kmer positions on negative node to equivalent bp positions?
//TODO: eliminate the need to use the SequenceMatch class
//    instead use an index returning the matches.
//    Possible optimisation: matches as const pointers to the index internal vectors rather than copying them
//    this may need to deal with FW/RC in a different way, but could be faster! (or would it just waste locality?)


class PerfectMatch{
public:
    PerfectMatch(sgNodeID_t _node=0,uint32_t _node_position=0,uint32_t _read_position=0,uint16_t _size=0): node(_node),node_position(_node_position),read_position(_read_position),size(_size){};
    //uint64_t read_id; //Not needed because we're using a vector per read.
    sgNodeID_t node;
    uint32_t node_position; //position is the start of the match on positive-node-coordinates
    uint32_t read_position; //position is the start of the match on positive-node-coordinates
    uint16_t size;

};

class LongReadsRecruiter {
public:
    LongReadsRecruiter(const SequenceDistanceGraph &sdg, const LongReadsDatastore &datastore,uint8_t k=25, uint16_t f=50);
    void recruit_reads(uint16_t seed_size,uint16_t seed_count,uint64_t first_read=1,uint64_t last_read=0);
    void reset_recruitment();
    //TODO: dump, load
    NKmerIndex nkindex;
    const SequenceDistanceGraph & sdg;
    const LongReadsDatastore &datastore;
    uint8_t k;
    uint16_t f;
    std::vector<std::vector<PerfectMatch>> read_perfect_matches;
    std::vector<std::vector<uint64_t>> node_reads;
};


#endif //SDG_LONGREADSRECRUITER_HPP
