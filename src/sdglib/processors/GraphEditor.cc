//
// Created by Bernardo Clavijo (EI) on 25/06/2018.
//

#include <functional>
#include "GraphEditor.hpp"

GraphEditorNodeDeletion::GraphEditorNodeDeletion(sgNodeID_t node, const std::vector<sgNodeID_t> &connected_ends){
    input_nodes.emplace_back(node);
    consumed_nodes=input_nodes;
    consumed_ends=connected_ends;
}

GraphEditorLinkDeletion::GraphEditorLinkDeletion(sgNodeID_t src, sgNodeID_t dest) {
    input_ends.emplace_back(src);
    input_ends.emplace_back(dest);
    consumed_ends=input_ends;
}

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

GraphEditorPathDetachment::GraphEditorPathDetachment(std::vector<sgNodeID_t> nodes,bool full) {
    full_detach=full;
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

bool GraphEditor::queue_node_deletion(sgNodeID_t node) {
    std::vector<sgNodeID_t> connected_ends(ws.sdg.links[llabs(node)].size());
    for (auto l:ws.sdg.links[llabs(node)]) connected_ends.emplace_back(l.dest);
    GraphEditorNodeDeletion op(node,connected_ends);
    //if (not queue_allows(op)) return false;
    //TODO: do we validate that the expansion is valid?
    queue_mark_inputs(op);
    op.index=next_op++;
    node_deletion_queue.emplace_back(op);
    return true;
}

bool GraphEditor::queue_link_deletion(sgNodeID_t src, sgNodeID_t dest) {
    GraphEditorLinkDeletion op(src,dest);
    //if (not queue_allows(op)) return false;
    //TODO: do we validate that the expansion is valid?
    queue_mark_inputs(op);
    op.index=next_op++;
    link_deletion_queue.emplace_back(op);
    return true;
}

bool GraphEditor::queue_node_expansion(sgNodeID_t node, std::vector<std::pair<sgNodeID_t, sgNodeID_t>> links_through) {
    GraphEditorNodeExpansion op(node,links_through);
    if (not queue_allows(op)) return false;
    //TODO: do we validate that the expansion is valid?
    queue_mark_inputs(op);
    op.index=next_op++;
    node_expansion_queue.emplace_back(op);
    return true;
}

bool GraphEditor::queue_path_detachment(std::vector<sgNodeID_t> nodes, bool full) {
    //Check path is valid
    if (nodes.size()<2) return true;
    SequenceDistanceGraphPath p(ws.sdg,nodes);
    if (p.is_unitig()) return true;
    if (!p.is_valid()) return false;
    GraphEditorPathDetachment op(nodes,full);
    if (not queue_allows(op)) return false;
    //TODO: do we validate that the expansion is valid?
    queue_mark_inputs(op);
    op.index=next_op++;
    path_detachment_queue.emplace_back(op);
    return true;

    //Mark first and last node INTERNAL ends as used
    //for (auto )

}

void GraphEditor::apply_all() {
    auto next_node_deletion=node_deletion_queue.begin();
    auto next_link_deletion=link_deletion_queue.begin();
    auto next_expansion=node_expansion_queue.begin();
    auto next_detachment=path_detachment_queue.begin();
    for (uint64_t i=i;i<next_op;++i) {
        if (next_node_deletion!=node_deletion_queue.end() and next_node_deletion->index==i) {
            auto &op=*next_node_deletion;
            ws.sdg.remove_node(llabs(op.input_nodes[0]));
            ++next_node_deletion;
        }
        if (next_link_deletion!=link_deletion_queue.end() and next_link_deletion->index==i) {
            auto &op=*next_link_deletion;
            ws.sdg.remove_link(op.input_ends[0],op.input_ends[1]);
            ++next_link_deletion;
        }
        if (next_expansion!=node_expansion_queue.end() and next_expansion->index==i) {
            auto &op=*next_expansion;
            for (auto li = 0; li < op.input_ends.size() / 2; ++li) {
                //create a copy of the node;
                auto new_node = ws.sdg.add_node(ws.sdg.get_node_sequence(op.input_nodes[0]));
                //connect to the sides
                auto source_dist = ws.sdg.get_link(op.input_ends[2 * li], op.input_nodes[0]).dist;
                ws.sdg.add_link(op.input_ends[2 * li], new_node, source_dist);
                auto dest_dist = ws.sdg.get_link(-op.input_nodes[0], op.input_ends[2 * li + 1]).dist;
                ws.sdg.add_link(-new_node, op.input_ends[2 * li + 1], dest_dist);
            }
            ws.sdg.remove_node(op.input_nodes[0]);
            ++next_expansion;
        }
        if (next_detachment!=path_detachment_queue.end() and next_detachment->index==i) {
            auto &op=*next_detachment;
            //special case: path with just ends
            if (op.input_nodes.empty()) {
                if ( op.full_detach ) {
                    for (auto l:ws.sdg.get_fw_links(-op.input_ends[0]))
                        if (l.dest != op.input_ends[1]) ws.sdg.remove_link(l.source, l.dest);
                    for (auto l:ws.sdg.get_bw_links(op.input_ends[1]))
                        if (l.dest != op.input_ends[0]) ws.sdg.remove_link(l.source, l.dest);

                }
            } else {
                //Create a unitig with the internal path sequence
                auto new_node = ws.sdg.add_node(SequenceDistanceGraphPath(ws.sdg, op.input_nodes).sequence());
                auto first_dist = ws.sdg.get_link(op.input_ends[0], op.input_nodes.front()).dist;
                auto last_dist = ws.sdg.get_link(-op.input_nodes.back(), op.input_ends[1]).dist;
                if ( op.full_detach ) {
                    for (auto l:ws.sdg.get_fw_links(-op.input_ends[0])) ws.sdg.remove_link(l.source, l.dest);
                    for (auto l:ws.sdg.get_bw_links(op.input_ends[1])) ws.sdg.remove_link(l.source, l.dest);
                }
                else {
                    ws.sdg.remove_link(op.input_ends[0], op.input_nodes.front());
                    ws.sdg.remove_link(-op.input_nodes.back(), op.input_ends[1]);
                }
                ws.sdg.add_link(op.input_ends[0], new_node, first_dist);
                ws.sdg.add_link(-new_node, op.input_ends[1], last_dist);
            }
            ++next_detachment;
        }

    }
    node_expansion_queue.clear();
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
        if (total>max_total) continue;
        else to_remove.insert(to_remove.end(),c.begin(),c.end());
    }
    uint64_t tbp=0;
    for (auto n:to_remove) tbp+=ws.sdg.nodes[llabs(n)].sequence.size();
    std::cout<<"There are "<<to_remove.size()<<" nodes and "<<tbp<<"bp in small unconnected components"<<std::endl;
    for (auto n:to_remove) ws.sdg.remove_node(llabs(n));
}