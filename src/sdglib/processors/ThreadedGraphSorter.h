//
// Created by Gonzalo Garcia (EI) on 2020-11-10.
//

#ifndef SDG_THREADEDGRAPHSORTER_H
#define SDG_THREADEDGRAPHSORTER_H


#include <sdglib/views/NodeView.hpp>
std::array<uint64_t,3> assess_node_happiness(sgNodeID_t nid, const std::unordered_map<sgNodeID_t , uint32_t> &order, const DistanceGraph& trg_nt);
std::map<sgNodeID_t , int64_t > sort_cc(const DistanceGraph& dg, std::unordered_set<sgNodeID_t> cc);
bool pop_node(DistanceGraph& dg, sgNodeID_t node_id, uint64_t read);
void pop_node_from_all(DistanceGraph& dg, sgNodeID_t nid);

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

class TheGreedySorter {
public:
    TheGreedySorter(const DistanceGraph& _trg_nt, sgNodeID_t founding_node=0);
    void update_read_nodes_in_order();
    std::vector<uint64_t> node_belonging_scores(uint64_t nid);
    std::vector<uint64_t > rids_from_node(NodeView nv);
    uint64_t shared_reads(NodeView nv1, NodeView nv2);
    bool pop_node(sgNodeID_t node_id, uint64_t read);
    bool remove_node(sgNodeID_t node_id);
    void start_from_read(uint64_t rid, int min_confirmation);
    std::pair<int, int> evaluate_read(uint64_t rid, bool print_pos);
    std::vector<int64_t> thread_node_positions(int64_t rid);
    std::vector<int32_t > evaluate_read_nodeorder(uint64_t rid, bool print_pos);
    void add_read(uint64_t rid, int min_confirmation=2);
//    void TheGreedySorter::extend_solution(int min_support=2, int min_shared=10, int min_new=10);

    void write_connected_nodes_graph(std::string filename);

    std::map<sgNodeID_t , int64_t > sort_graph();

    const DistanceGraph& trg_nt;
    DistanceGraph dg;

    std::vector<sgNodeID_t > all_nodes;
    std::unordered_set<uint64_t > all_reads;
    std::set<sgNodeID_t > used_nodes;
    std::vector<sgNodeID_t > used_reads;
    std::map<uint64_t, std::vector<sgNodeID_t >> read_ends;

private:
    bool order_is_valid=false;
    std::map<sgNodeID_t , int64_t > order;
    std::map<uint64_t,uint64_t> read_nodes_in_order;
};

#endif //SDG_THREADEDGRAPHSORTER_H
