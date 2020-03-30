//
// Created by Bernardo Clavijo (EI) on 25/06/2018.
//

#include <functional>
#include "GraphEditor.hpp"

GraphEditorNodeExpansion::GraphEditorNodeExpansion(sgNodeID_t node,
                                                   std::vector<std::pair<sgNodeID_t, sgNodeID_t>> links_through) {

    input_nodes.emplace_back(node);
    for (auto &l:links_through) {
        input_ends.emplace_back(l.first);
        input_ends.emplace_back(l.second);
    }
    consumed_nodes=input_nodes;
    consumed_ends=input_ends;
}

GraphEditorPathDetachment::GraphEditorPathDetachment(std::vector<sgNodeID_t> nodes) {
    //Check first and last have their INTERNAL ends available
    input_ends.emplace_back(-nodes.front());
    input_ends.emplace_back(nodes.back());
    consumed_ends=input_ends;
    //Check every node in the path is FULLY AVAILABLE
    for (auto i=1;i<nodes.size()-1;++i)
        input_nodes.emplace_back(nodes[i]);
}

bool GraphEditor::queue_allows(GraphEditorOperation op) {
    for (auto &n:op.input_nodes) {
        if (queued_nodes[llabs(n)] or queued_plus_ends[llabs(n)] or queued_minus_ends[llabs(n)]) return false;
    }
    for (auto &n:op.input_ends) {
        if (queued_nodes[llabs(n)]) return false;
        if (n < 0 and queued_minus_ends[-n]) return false;
        if (n > 0 and queued_plus_ends[n]) return false;
    }
    return true;
}

void GraphEditor::queue_mark_inputs(GraphEditorOperation op) {
    for (auto &n:op.consumed_nodes) {
        queued_nodes[llabs(n)]=true;
        queued_plus_ends[llabs(n)]=true;
        queued_minus_ends[llabs(n)]=true;
    }
    for (auto &e:op.consumed_ends) {
        if (e < 0) queued_minus_ends[-e]=true;
        if (e > 0) queued_plus_ends[e]=true;
    }
}

bool GraphEditor::queue_node_expansion(sgNodeID_t node, std::vector<std::pair<sgNodeID_t, sgNodeID_t>> links_through) {
    GraphEditorNodeExpansion op(node,links_through);
    if (not queue_allows(op)) return false;
    //TODO: do we validate that the expansion is valid?
    queue_mark_inputs(op);
    node_expansion_queue.emplace_back(op);
    return true;
}

bool GraphEditor::queue_path_detachment(std::vector<sgNodeID_t> nodes) {
    //Check path is valid
    if (nodes.size()<2) return true;
    SequenceDistanceGraphPath p(ws.sdg,nodes);
    if (p.is_unitig()) return true;
    if (!p.is_valid()) return false;
    GraphEditorPathDetachment op(nodes);
    if (not queue_allows(op)) return false;
    //TODO: do we validate that the expansion is valid?
    queue_mark_inputs(op);
    path_detachment_queue.emplace_back(op);
    return true;

    //Mark first and last node INTERNAL ends as used
    //for (auto )

}

void GraphEditor::apply_all(bool remove_small_components_total_bp) {
    for (auto &op:node_expansion_queue){
        for (auto li=0;li<op.input_ends.size()/2;++li) {
            //create a copy of the node;
            auto new_node=ws.sdg.add_node(ws.sdg.get_node_sequence(op.input_nodes[0]));
            //connect to the sides
            auto source_dist=ws.sdg.get_link(op.input_ends[2*li],op.input_nodes[0]).dist;
            ws.sdg.add_link(op.input_ends[2*li],new_node,source_dist);
            auto dest_dist=ws.sdg.get_link(-op.input_nodes[0],op.input_ends[2*li+1]).dist;
            ws.sdg.add_link(-new_node,op.input_ends[2*li+1],dest_dist);
        }
        ws.sdg.remove_node(op.input_nodes[0]);
    }
    node_expansion_queue.clear();
    for (auto &op:path_detachment_queue){
        std::cout<<"applying "<<-op.input_ends[0]<<" -> (";
        for (auto &n:op.input_nodes) std::cout<<" "<<n;
        std::cout<<" ) -> "<<op.input_ends[1]<<std::endl;
        //special case: path with just ends
        if (op.input_nodes.empty()){
            for (auto l:ws.sdg.get_fw_links(-op.input_ends[0])){
                if (l.dest!=op.input_ends[1])
                    ws.sdg.remove_link(l.source,l.dest);
            }
        }
        else {
            //Create a unitig with the internal path sequence
            auto new_node = ws.sdg.add_node(SequenceDistanceGraphPath(ws.sdg, op.input_nodes).sequence());
            std::cout<<" new node="<<new_node<<std::endl;
            auto first_dist = ws.sdg.get_link(op.input_ends[0], op.input_nodes.front()).dist;
            auto last_dist = ws.sdg.get_link(-op.input_nodes.back(), op.input_ends[1]).dist;
            //for (auto l:ws.sdg.get_fw_links(-op.input_ends[0])) ws.sdg.remove_link(l.source, l.dest);
            ws.sdg.remove_link(op.input_ends[0],op.input_nodes.front());
            //for (auto l:ws.sdg.get_bw_links(op.input_ends[1])) ws.sdg.remove_link(l.source, l.dest);
            ws.sdg.remove_link(-op.input_nodes.back(),op.input_ends[1]);
            ws.sdg.add_link(op.input_ends[0], new_node, first_dist);
            ws.sdg.add_link(-new_node, op.input_ends[1], last_dist);
        }
    }
    path_detachment_queue.clear();
}


void GraphEditor::remove_small_components(int max_nodes, int max_size, int max_total) {
    std::vector<sgNodeID_t> to_remove;
    for (auto c:ws.sdg.connected_components(0,0,0)){
        if (c.size()>max_nodes) continue;
        uint64_t total=0;
        for (auto n:c) {
            if (ws.sdg.nodes[llabs(n)].sequence.size()>max_size) total+=max_total;
            total+=ws.sdg.nodes[llabs(n)].sequence.size();
        }
        if (total>max_size) continue;
        else to_remove.insert(to_remove.end(),c.begin(),c.end());
    }
    uint64_t tbp=0;
    for (auto n:to_remove) tbp+=ws.sdg.nodes[llabs(n)].sequence.size();
    std::cout<<"There are "<<to_remove.size()<<" nodes and "<<tbp<<"bp in small unconnected components"<<std::endl;
    for (auto n:to_remove) ws.sdg.remove_node(llabs(n));
}