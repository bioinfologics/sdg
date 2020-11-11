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

    bool pop_node(sgNodeID_t node_id, uint64_t read);

    bool multipop_node(sgNodeID_t node_id, uint64_t read);

    // Not implemented bc not used in the pipeline
//    DistanceGraph make_summary_connection_graph(DistanceGraph& mdg, int links_to_connect=3){};

    // TODO: this is not filtering the nodes, for some reason refuses to take the connected only parameter in get_all_nodeviews
    void write_connected_nodes_graph(std::string filename);

    std::pair<std::map<sgNodeID_t , int64_t >, std::vector<sgNodeID_t >> sort_cc(std::vector<sgNodeID_t> cc);

    const DistanceGraph & dg;
};


#endif //SDG_THREADEDGRAPHSORTER_H
