//
// Created by Bernardo Clavijo (EI) on 16/01/2020.
//

#ifndef SDG_LONGREADSRECRUITER_HPP
#define SDG_LONGREADSRECRUITER_HPP


#include <sdglib/graph/SequenceDistanceGraph.hpp>
#include <sdglib/datastores/LongReadsDatastore.hpp>

//TODO: OMP
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

class NodePosition{
public:
    NodePosition(sgNodeID_t _node,int32_t _start,int32_t _end): node(_node),start(_start),end(_end){};
    //uint64_t read_id; //Not needed because we're using a vector per read.
    bool operator<(const struct NodePosition &other) const{
        return std::tie(start,end,node)<std::tie(other.start,other.end,other.node);
    }
    sgNodeID_t node;
    int32_t start; //we need both start and end due to read indels
    int32_t end; //we need both start and end due to read indels
};

class LongReadsRecruiter {
public:
    LongReadsRecruiter(SequenceDistanceGraph &sdg, const LongReadsDatastore &datastore,uint8_t k=25, uint16_t f=50);
    void dump(std::string filename);
    void load(std::string filename);
    void perfect_mappings(uint16_t seed_size,uint64_t first_read=1,uint64_t last_read=0);
    void recruit_reads(uint16_t seed_size,uint16_t seed_count,uint64_t first_read=1,uint64_t last_read=0);
    void reset_recruitment();

    std::vector<NodePosition> endmatches_to_positions(uint64_t rid,int32_t end_size, uint16_t matches);
    void thread_reads(uint32_t end_size, uint16_t matches); //uses endmatches_to_positions
    DistanceGraph dg_from_threads(bool multi_link=false);

    SequenceDistanceGraph & sdg;
    const LongReadsDatastore &datastore;
    uint8_t k;
    uint16_t f;
    std::vector<std::vector<PerfectMatch>> read_perfect_matches;
    std::vector<std::vector<NodePosition>> read_threads;
    std::vector<std::vector<uint64_t>> node_reads;
};


#endif //SDG_LONGREADSRECRUITER_HPP
