//
// Created by Bernardo Clavijo (EI) on 16/01/2020.
//

#pragma once


#include <sdglib/graph/SequenceDistanceGraph.hpp>
#include <sdglib/datastores/LongReadsDatastore.hpp>
class ReadThreadsGraph;//fw declaration

//TODO: OMP
//TODO: translate kmer positions on negative node to equivalent bp positions?
//TODO: eliminate the need to use the SequenceMatch class
//    instead use an index returning the matches.
//    Possible optimisation: matches as const pointers to the index internal vectors rather than copying them
//    this may need to deal with FW/RC in a different way, but could be faster! (or would it just waste locality?)


/** @brief Describes a perfect match between 2 sequences (usually a node and a read)
 *
 * within a perfect match bot sequences are EXACTLY the same, no exceptions!.
 *
 */
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

/**
 * Class to store a node position in a read.
 * Used to describe mapping results
 */
class NodePosition{
public:
    NodePosition(sgNodeID_t _node=0,int32_t _start=0,int32_t _end=0): node(_node),start(_start),end(_end){};
    //uint64_t read_id; //Not needed because we're using a vector per read.
    bool operator<(const struct NodePosition &other) const{
        return std::tie(start,end,node)<std::tie(other.start,other.end,other.node);
    }
    sgNodeID_t node;
    int32_t start; //we need both start and end due to read indels
    int32_t end; //we need both start and end due to read indels
};

/**
 * Class to thread long reads/sequences to graphs
 */
class LongReadsRecruiter {
public:
    /** @brief LongReadsRecruiter initialization
     *
     * Sert internal variables and clean any previous recruitment
     * @param sdg Base graph
     * @param datastore long reads datastore handle to map
     * @param k mapping k
     * @param f max frequency for the MkmerIndex index
     */
    LongReadsRecruiter(SequenceDistanceGraph &sdg, const LongReadsDatastore &datastore,uint8_t k=25, uint16_t f=50);

    /** @brief dumps reads perfect matches mappings to a file for persistence
     *  To be restored using load
     *
     * @param filename output filename
     */
    void dump(std::string filename);

    /** @brief dumps read threads to a file for persistence
     *  To be restored using  load_threads
     *
     * @param filename output filename
     */
    void dump_threads(std::string filename);
    void load(std::string filename);
    void load_threads(std::string filename);

    /** @brief map reads by perfect matching stretches
     *
     *  Fills read_perfect_matches ( std::vector<std::vector<PerfectMatch>> ) with matches as they appear in the read.
     *
     *  Matches don't span multiple nodes, Nodes need to be of at least seed_size to generate matches.
     *
     * @param seed_size Minimum kmer stretch for a match to be considered
     * @param first_read position in the datastore of the first read to map (chunk initial position)
     * @param last_read positino in the datastore the last read to map (chunk final position)
     */
    void perfect_mappings(uint16_t seed_size,uint64_t first_read=1,uint64_t last_read=0);
    std::vector<PerfectMatch> reverse_perfect_matches(const std::vector<PerfectMatch> &matches, uint64_t rsize=0) const;

    /**
     * Maps the read by producing PerfectMatches and storing the results in read_perfect_matches.
     * This function is capable of continuing a match over a link (across overlap), thus small nodes are covered by a
     * match and can produce usable mappings (the matches are stored independently by node anyway)
     *
     * Fills read_perfect_matches ( std::vector<std::vector<PerfectMatch>> ) with perfect matches per read
     *
     * @param seed_size Minimum match size
     * @param first_read position in the datastore of the first read to map (chunk initial position)
     * @param last_read positino in the datastore the last read to map (chunk final position)
     */
    void map(uint16_t seed_size);

    void anchormap(uint16_t seed_size, std::string kcname, std::string countname, uint32_t fmin, uint32_t fmax, uint32_t graph_fmin=1, uint32_t graph_fmax=1);

    void map_to_index(uint16_t seed_size, NKmerIndex &nki);

    /** @brief recruit the reads in nodes.
     * Fills the node_reads ( std::vector<std::vector<int64_t>> ) vector.
     *
     * Takes all matches and assigns each match to the corresponding node,
     *
     * the node_reads[node_id] position stores all reads where node_id was mapped to.
     * to get the mapped reads to the node.
     * Contrary to recruit_threads() this function uses read_perfect_matches instead of the read_threads
     *
     * read_perfect_matches can be filled using map()
     *
     * @param seed_size Minimum seed size for a match to be recruited
     * @param seed_count Number of matches for a read to be recruited
     * @param first_read position in the datastore of the first read to map (chunk initial position)
     * @param last_read positino in the datastore the last read to map (chunk final position)
     */
    void recruit_reads(uint16_t seed_size,uint16_t seed_count,int64_t first_read=1,int64_t last_read=0);

    /** @brief recruits the threads in the nodes.
     * Fills the node_threads ( std::vector<std::vector<int64_t>> ) vector.
     *
     * For each node all threads that go over are recruited.
     * The node_threads[node_id] position stores all threads where node_id is included.
     * Contrary to recruit_reads() this function uses read_threads instead of the read_perfect_matches
     *
     * read_threads can be filled using thread_reads(), simple_thread_reads(), thread_and_pop()
     */
    void recruit_threads();

    /** @brief Resets all recruitment
     * Erase and resize read_perfect_matches and read_threads
     */
    void reset_recruitment();

    /** @brief get a list of nodes fw from node in read_id
     *
     * If a negative node is used the function will return the BW path (the FW of the negative node)
     *
     * @param read_id read to search the fw path in
     * @param node node to search the fw path from (neg for BW path)
     * @return list of FW nodes from node in read_id
     */
    std::vector<sgNodeID_t> path_fw(seqID_t read_id, sgNodeID_t node) const;

    /** @brief get all fw paths from a node
     *
     * Similar to path_fw() but returns fw paths for all available reads in the node_paths vector
     *
     * @param node node to get the paths fw from
     * @return all paths fw from the node
     */
    std::vector<std::vector<sgNodeID_t> > all_paths_fw(sgNodeID_t node) const;


    std::pair<int,int> find_next_valid_block(const std::vector<PerfectMatch> & matches, int start, int min_count);
    /** @brief Threads the nodes by plain aggregation of the matches to a node
     *
     * Matches don't need to fulfill any other requirement to appear in the thread other thatn to have a match
     * Fills read_threads with NodePositions. read_threads[i] stores the NodePositions of the thread created by read i.
     *
     */
    void simple_thread_reads(int min_count=1);

    /** @brief Create a ReadThreadsGraph object using the threads information
     * see ReadThreadsGraph class
     * @param remove_duplicated remove duplicated links
     * @param min_thread_nodes Min length of the thread in nodes for the thread to be considered
     * @return ReadThreadsGraph object with all the threads in read_threads added
     */
    ReadThreadsGraph rtg_from_threads(bool remove_duplicated=false, int min_thread_nodes=1);

    SequenceDistanceGraph & sdg;
    const LongReadsDatastore &datastore;
    uint8_t k;
    uint16_t f;
    std::vector<std::vector<PerfectMatch>> read_perfect_matches;
    std::vector<std::vector<NodePosition>> read_threads;
    std::vector<std::vector<int64_t>> node_reads; //TODO: sgNodeID_t represented as int64_t, shoud not couse problems but WHY?
    std::vector<std::vector<int64_t>> node_threads;
    std::vector<std::vector<int32_t>> node_paths;
    std::vector<std::vector<sgNodeID_t>> read_paths;
};



