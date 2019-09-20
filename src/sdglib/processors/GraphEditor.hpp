//
// Created by Bernardo Clavijo (EI) on 25/06/2018.
//

#pragma once

#include <sdglib/workspace/WorkSpace.hpp>
class GraphEditorOperation {
public:
    std::vector<sgNodeID_t> input_nodes; //nodes need to be fully untouched by queue
    std::vector<Link> input_links; //ends in links need to be fully untouched by queue
//    std::vector<sgNodeID_t> consumed_nodes;
//    std::vector<Link> consumed_links;
//    std::vector<sgNodeID_t> generated_nodes;
//    std::vector<Link> generated_links;
};

/**
 * An operation expanding a node into as many copies as links-through.
 * For every linkf_through=(s,d) a copy nc is create and links (s,nc) and (-nc,d) are created.
 */
class GraphEditorNodeExpansion : public GraphEditorOperation {
public:
    GraphEditorNodeExpansion(sgNodeID_t node,std::vector<std::pair<sgNodeID_t,sgNodeID_t>> links_through);
//    std::vector<std::pair<sgNodeID_t,sgNodeID_t>> links_through;
};

class GraphEditor {
public:
    GraphEditor (WorkSpace &_ws):ws(_ws){
        queued_nodes.resize(ws.sdg.nodes.size());
        queued_plus_ends.resize(ws.sdg.nodes.size());
        queued_minus_ends.resize(ws.sdg.nodes.size());
    };

    /**
     * Puts a node expansion in the queue.
     * @param node
     * @param links_through
     * @return true if operation added, false if conflict with queue detected
     */
    bool queue_node_expansion(sgNodeID_t node,std::vector<std::pair<sgNodeID_t,sgNodeID_t>> links_through);

    /**
     * Checks if an operation is allowed given the current queue.
     * @return true if operation allowed
     */
    bool queue_allows(GraphEditorOperation);

    /**
    * Marks inputs as used in the queue.
    * @return true if operation allowed
    */
    bool queue_mark_inputs(GraphEditorOperation);

    /**
     * applies the queued operations
     */
    void apply_all();

    /******** OLD DEPRECATED FUNCTIONS *******/
    bool detach_path (SequenceDistanceGraphPath p,bool consume_tips=false); //returns true on success
    int patch_between (sgNodeID_t from, sgNodeID_t to, std::string);
    SequenceDistanceGraphPath find_longest_path_from(sgNodeID_t node, std::string seq);
    void remove_small_components(int max_nodes, int max_size, int max_total);

    /**
     * @brief expands a node creating as many copies as needed, then, distributes input and output links as per bw and fw
     * @param nodeID
     * @param bw
     * @param fw
     */
    void expand_node(sgNodeID_t nodeID, std::vector<std::vector<sgNodeID_t>> bw, std::vector<std::vector<sgNodeID_t>> fw);

    void expand_path(const SequenceDistanceGraphPath &p);

    void join_path(SequenceDistanceGraphPath p,bool consume_nodes=false);

    WorkSpace & ws;
    std::set<sgNodeID_t> edited_nodes;
    std::set<sgNodeID_t> new_nodes;

    //Node and node end status are saved as bool vectors for fast lookup
    std::vector<bool> queued_nodes;
    std::vector<bool> queued_plus_ends;
    std::vector<bool> queued_minus_ends;

    //Queues are by operation type, order does not matter as all operations are compatible
    std::vector<GraphEditorNodeExpansion> node_expansion_queue;

};

