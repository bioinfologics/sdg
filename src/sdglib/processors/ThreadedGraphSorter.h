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

    // De la otra celda
    std::vector<uint64_t > rids_from_node(NodeView nv);
//    shared_reads
//    pop_node
//    multipop_node
//    make_summary_connection_graph
//    write_connected_nodes_graph
//    sort_cc
    const DistanceGraph & dg;
};


#endif //SDG_THREADEDGRAPHSORTER_H
