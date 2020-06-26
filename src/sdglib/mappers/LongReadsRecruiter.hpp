//
// Created by Bernardo Clavijo (EI) on 16/01/2020.
//

#pragma once


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
    PerfectMatch(sgNodeID_t _node=0,int32_t _node_position=0,int32_t _read_position=0,uint16_t _size=0): node(_node),node_position(_node_position),read_position(_read_position),size(_size){};
    //uint64_t read_id; //Not needed because we're using a vector per read.
    bool operator<(const struct PerfectMatch &other) const{
        return std::tie(read_position,size,node,node_position)<std::tie(other.read_position,other.size,other.node,other.node_position);
    }
    sgNodeID_t node;
    int32_t node_position; //position is the start of the match on positive-node-coordinates
    int32_t read_position; //position is the start of the match on positive-node-coordinates
    uint16_t size;

};

class PerfectMatchesFilter {
public:
    PerfectMatchesFilter( WorkSpace &_ws):ws(_ws){};
    std::vector<PerfectMatch> truncate_turnaround (const std::vector<PerfectMatch> &in) const;
    std::vector<PerfectMatch> matches_fw_from_node (sgNodeID_t node, const std::vector<PerfectMatch> &in) const;
    std::vector<PerfectMatch> clean_linear_groups(const std::vector<PerfectMatch> &in, int group_size=5,int small_node_size=500) const;
    std::vector<PerfectMatch> merge_and_sort(const std::vector<std::vector<PerfectMatch>> &in) const;
    WorkSpace & ws;
};
class LongReadsRecruiter;

class PerfectMatchesMergeSorter{
public:
    PerfectMatchesMergeSorter( WorkSpace &_ws):ws(_ws){};

    void init_from_node(sgNodeID_t n,const LongReadsRecruiter &lrr, int min_reads=3, int group_size=5, int small_node_size=500);
    void drop_conflictive_reads();
    void find_next_node(int d=2000, float candidate_percentaje=.6, float first_percentaje=.95, bool verbose=false);
    void advance_reads_to_node();
    void advance_reads_through_node();
    std::vector<std::vector<PerfectMatch>> read_matches;
    std::vector<int32_t> read_next_match; //-1 means read has been dropped/finished
    std::vector<int32_t> read_last_hit_position;
    std::vector<int32_t> read_dropped_position; //-1 means read has not been dropped
    std::vector<PerfectMatch> out;
    sgNodeID_t next_node;
    WorkSpace & ws;
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
    std::vector<PerfectMatch> reverse_perfect_matches(const std::vector<PerfectMatch> &matches, uint64_t rsize=0);
    void map(uint16_t seed_size,uint64_t first_read=1,uint64_t last_read=0);
    void recruit_reads(uint16_t seed_size,uint16_t seed_count,int64_t first_read=1,int64_t last_read=0);
    void recruit_threads();
    void reset_recruitment();
    void thread_and_pop();
    std::vector<sgNodeID_t> path_fw(seqID_t read_id, sgNodeID_t node) const;
    std::vector<std::vector<sgNodeID_t> > all_paths_fw(sgNodeID_t node) const;

    std::vector<NodePosition> endmatches_to_positions(uint64_t rid,int32_t end_size, uint16_t matches);
    void thread_reads(uint32_t end_size, uint16_t matches); //uses endmatches_to_positions
    DistanceGraph dg_from_threads(bool multi_link=false);

    SequenceDistanceGraph & sdg;
    const LongReadsDatastore &datastore;
    uint8_t k;
    uint16_t f;
    std::vector<std::vector<PerfectMatch>> read_perfect_matches;
    std::vector<std::vector<NodePosition>> read_threads;
    std::vector<std::vector<int64_t>> node_reads;
    std::vector<std::vector<int64_t>> node_threads;
    std::vector<std::vector<int32_t>> node_paths;
    std::vector<std::vector<sgNodeID_t>> read_paths;
};


