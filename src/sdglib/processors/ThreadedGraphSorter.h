//
// Created by Gonzalo Garcia (EI) on 2020-11-10.
//

#ifndef SDG_THREADEDGRAPHSORTER_H
#define SDG_THREADEDGRAPHSORTER_H


#include <sdglib/views/NodeView.hpp>
#include <sdglib/mappers/LongReadsRecruiter.hpp>
std::array<uint64_t,3> assess_node_happiness(sgNodeID_t nid, const std::unordered_map<sgNodeID_t , uint32_t> &order, const DistanceGraph& trg_nt);
std::unordered_map<sgNodeID_t , int64_t > sort_cc(const DistanceGraph& dg, std::unordered_set<sgNodeID_t> cc);
bool pop_node(DistanceGraph& dg, sgNodeID_t node_id, uint64_t read);
void pop_node_from_all(DistanceGraph& dg, sgNodeID_t nid);


//This takes a thread and pops all unhappy/disconnected nodes from it, returns an empty thread if theres too many of them
std::vector<NodePosition> make_thread_happy(const std::vector<NodePosition> &thread,const DistanceGraph & trg, int max_unhappy=1, float disconnection_rate=.3);
void make_all_threads_happy(LongReadsRecruiter & lrr, DistanceGraph &trg, int max_unhappy=1, float disconnection_rate=.3);

//class ThreadedGraphSorter {
//public:
//    // TODO read ids (rids) are uint64_t type always, should we define a rid type?? YES
//    ThreadedGraphSorter(const DistanceGraph & _dg){};
//
//    /*
//     * Returns a list of reads linking to or from the spec nodeview
//     * TODO: Change nv for node maybe!?
//     * */
////    std::vector<uint64_t > rids_from_node(NodeView nv);
////    uint64_t shared_reads(NodeView nv1, NodeView nv2);
//
////    bool pop_node(sgNodeID_t node_id, uint64_t read);
////    bool multipop_node(sgNodeID_t node_id, uint64_t read);
//
//    // Not implemented bc not used in the pipeline
////    DistanceGraph make_summary_connection_graph(DistanceGraph& mdg, int links_to_connect=3){};
//
//
//
//    const DistanceGraph & dg;
//};

enum HappyPosition{Nowhere,FrontFW,FrontBW,MiddleFW,MiddleBW,BackFW,BackBW};
class HappyInsertionSorter;//FW declaration
/**
 * this keeps a record of a Node's pre/post conditions and computes if the node is happy to be added to the order
 */
class NodeAdjacencies {
public:
    NodeAdjacencies (){};
    void mark_used(const sgNodeID_t & nid);
    void mark_unused(const sgNodeID_t & nid);
    HappyPosition happy_to_add(float used_perc); //TODO::expand to return a HappyPosition
    std::pair<int64_t,int64_t> find_happy_place(const HappyInsertionSorter& sorter); //finds the last of the prevs and first of the nexts, signed for orientation
    std::unordered_set<sgNodeID_t> prevs,nexts;
    uint64_t fw_prevs=0,bw_prevs=0,fw_nexts=0,bw_nexts=0;
};
/**
 * insertion sorter: needs to be founded by a relatively good set of ordered nodes (i.e. the local order of a node)
 * keeps a record of previous/next nodes to every node: if most previous/next nodes to a node are in the order, add the node to it
 */
class HappyInsertionSorter {
public:
    HappyInsertionSorter(ReadThreadsGraph& _rtg):rtg(_rtg){};
    void compute_adjacencies(int min_links=2, int radius=20);
    void add_node(sgNodeID_t nid); //figures out orientation, then adds with sign, then runs mark_used on all this node's NodeAdjacencies (only those will be affected), removes this node from candidates and adds its adjacencies
    NodeAdjacencies & get_node_adjacencies(sgNodeID_t nid);
    std::vector<sgNodeID_t> nodes_to_add(float used_perc=.9); //finds nodes to add, sorts them by inter-dependency (i.e. nodes coming first are not made happier by subsequent nodes).
    bool insert_node(sgNodeID_t nid,float used_perc=0.9, bool solve_floating_by_rtg=false);//uses the node's happy place to put it in the order, shifts other nodes around if needed. Returns false if the node is unhappy.
    //TODO: validate_order -> checks all nodes are still happy in their place.
    void reset_positions();
    int64_t get_node_position(sgNodeID_t nid) const; //negative means reverse, but order is abs!!!!
    void remove_node_from_everywhere(sgNodeID_t nid);

    ReadThreadsGraph& rtg;
    std::unordered_set<sgNodeID_t> candidates;
    std::unordered_map<sgNodeID_t,int64_t> node_positions;
    std::unordered_map<sgNodeID_t,NodeAdjacencies> adjacencies;
};

class TheGreedySorter {
public:
    TheGreedySorter(const DistanceGraph& _trg_nt, sgNodeID_t founding_node=0);
    void update_read_nodes_in_order();
    std::vector<uint64_t> node_belonging_scores(int64_t nid);
    std::vector<uint64_t > rids_from_node(NodeView nv);
    uint64_t shared_reads(NodeView nv1, NodeView nv2);
    bool pop_node(sgNodeID_t node_id, uint64_t read);
    bool remove_node(sgNodeID_t node_id);
    void start_from_read(uint64_t rid, int min_confirmation);
    std::pair<sgNodeID_t,sgNodeID_t> get_thread_ends(int64_t rid);
    std::pair<int, int> evaluate_read(uint64_t rid, bool print_pos);
    std::vector<int64_t> thread_nodes(int64_t rid);
    std::vector<int64_t> thread_node_positions(int64_t rid);
    std::vector<int32_t > evaluate_read_nodeorder(uint64_t rid, bool print_pos);
    void add_read(uint64_t rid, int min_confirmation=2);
//    void TheGreedySorter::extend_solution(int min_support=2, int min_shared=10, int min_new=10);

    void write_connected_nodes_graph(std::string filename);

    std::unordered_map<sgNodeID_t , int64_t > sort_graph();

    const DistanceGraph& trg_nt;
    DistanceGraph dg;

    std::vector<sgNodeID_t > all_nodes;
    std::unordered_set<uint64_t > all_reads;
    std::set<sgNodeID_t > used_nodes;
    std::vector<sgNodeID_t > used_reads;
    std::unordered_map<uint64_t, std::vector<sgNodeID_t >> read_ends;

private:
    bool order_is_valid=false;
    std::unordered_map<sgNodeID_t , int64_t > order;
    std::unordered_map<uint64_t,uint64_t> read_nodes_in_order;
};

#endif //SDG_THREADEDGRAPHSORTER_H
