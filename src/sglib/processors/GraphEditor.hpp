//
// Created by Bernardo Clavijo (EI) on 25/06/2018.
//

#ifndef BSG_GRAPHEDITOR_HPP
#define BSG_GRAPHEDITOR_HPP


#include <sglib/workspace/WorkSpace.hpp>

class GraphEditor {
public:
    GraphEditor (WorkSpace &_ws):ws(_ws){};

    bool detach_path (SequenceGraphPath p,bool consume_tips=false); //returns true on success
    int patch_between (sgNodeID_t from, sgNodeID_t to, std::string);
    SequenceGraphPath find_longest_path_from(sgNodeID_t node, std::string seq);
    void remove_small_components(int max_nodes, int max_size, int max_total);

    /**
     * @brief expands a node creating as many copies as needed, then, distributes input and output links as per bw and fw
     * @param nodeID
     * @param bw
     * @param fw
     */
    void expand_node(sgNodeID_t nodeID, std::vector<std::vector<sgNodeID_t>> bw, std::vector<std::vector<sgNodeID_t>> fw);

    void expand_path(const SequenceGraphPath &p);

    void join_path(SequenceGraphPath p,bool consume_nodes=false);

    WorkSpace & ws;
    std::set<sgNodeID_t> edited_nodes;
    std::set<sgNodeID_t> new_nodes;

};


#endif //BSG_GRAPHEDITOR_HPP
