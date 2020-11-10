//
// Created by Gonzalo Garcia (EI) on 2020-11-10.
//

#ifndef SDG_THREADEDGRAPHSORTER_H
#define SDG_THREADEDGRAPHSORTER_H


#include <sdglib/views/NodeView.hpp>

class ThreadedGraphSorter {
public:
    // TODO read ids (rids) are uint64_t type always, should we define a rid type?? YES
    ThreadedGraphSorter(const DistanceGraph & _dg): dg(_dg){};
//    __init__
//    start_from_read
//    evaluate_read
//    evaluate_read_nodeorder
//    add_realted_nodes_and_reads
//    add_read
//    sort_graph

    /*
     * Returns a list of reads linking to or from the spec nodeview
     * TODO: Change nv for node maybe!?
     * */
    std::vector<uint64_t > rids_from_node(NodeView nv);
    uint64_t shared_reads(NodeView nv1, NodeView nv2);
//    pop_node
//    multipop_node
//    make_summary_connection_graph
//    write_connected_nodes_graph
//    sort_cc
    const DistanceGraph & dg;
};


#endif //SDG_THREADEDGRAPHSORTER_H
