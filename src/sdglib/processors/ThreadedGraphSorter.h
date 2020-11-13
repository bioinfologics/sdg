//
// Created by Gonzalo Garcia (EI) on 2020-11-10.
//

#ifndef SDG_THREADEDGRAPHSORTER_H
#define SDG_THREADEDGRAPHSORTER_H


#include <sdglib/views/NodeView.hpp>

std::map<sgNodeID_t , int64_t > sort_cc(const DistanceGraph& dg, std::unordered_set<sgNodeID_t> cc);

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
    TheGreedySorter(const DistanceGraph& _trg_nt);

    std::vector<uint64_t > rids_from_node(NodeView nv);
    uint64_t shared_reads(NodeView nv1, NodeView nv2);
    bool pop_node(sgNodeID_t node_id, uint64_t read);
    void start_from_read(uint64_t rid, int min_confirmation);
    std::pair<int, int> evaluate_read(uint64_t rid, bool print_pos);
    std::vector<int32_t > evaluate_read_nodeorder(uint64_t rid, bool print_pos);
    void add_read(uint64_t rid, int min_confirmation=2);
//    void TheGreedySorter::extend_solution(int min_support=2, int min_shared=10, int min_new=10);

    void write_connected_nodes_graph(std::string filename);

    std::map<sgNodeID_t , int64_t > sort_graph();

    const DistanceGraph& trg_nt;
    DistanceGraph dg;

    std::vector<sgNodeID_t > all_nodes;
    std::unordered_set<uint64_t > all_reads;
    std::vector<sgNodeID_t > used_nodes;
    std::vector<sgNodeID_t > used_reads;
    std::map<uint64_t, std::vector<sgNodeID_t >> read_ends;

private:
    bool order_is_valid=false;
    std::map<sgNodeID_t , int64_t > order;
};

#endif //SDG_THREADEDGRAPHSORTER_H
