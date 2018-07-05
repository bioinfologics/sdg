//
// Created by Bernardo Clavijo (EI) on 25/06/2018.
//

#include "GraphEditor.hpp"

bool GraphEditor::detach_path(SequenceGraphPath p, bool consume_tips) {
    SequenceGraph& sg(ws.getGraph());
    if (p.getNodes().size()==1) return true;
    //TODO: check the path is valid?
    if (p.getNodes().size()==2) {
        //std::cout<<"DETACHING path with only 2 nodes"<<std::endl;
        if (p.getNodes().size()==2){
            auto fwls=sg.get_fw_links(p.getNodes()[0]);
            for (auto l:fwls) if (l.dest!=p.getNodes()[1]) sg.remove_link(l.source,l.dest);
            auto bwls=sg.get_bw_links(p.getNodes()[1]);
            for (auto l:bwls) if (l.dest!=-p.getNodes()[0]) sg.remove_link(l.source,l.dest);
        }
        return true;
    }
    //check if the path is already detached.
    if (p.is_unitig()) return true;
    //check if the nodes in the path have already been modified
    for (auto n:p.getNodes()) {
        if (edited_nodes.count(llabs(n))>0) {
            std::cout<<"NOT DETACHING: path has already edited-out nodes"<<std::endl;
            return false;
        }
    }
    //TODO: create a new node with the sequence of the "middle" nodes
    auto pmid=p;
    pmid.getNodes().clear();
    for (auto i=1;i<p.getNodes().size()-1;++i) pmid.getNodes().push_back(p.getNodes()[i]);
    auto src=p.getNodes().front();
    auto dest=p.getNodes().back();
    if (!pmid.is_canonical()) {
        pmid.reverse();
        auto nsrc=-dest;
        dest=-src;
        src=nsrc;
    }
    sgNodeID_t new_node=sg.add_node(Node(pmid.get_sequence()));
    //TODO: migrate connection from the src node to the first middle node into the new node
    auto old_link=sg.get_link(-p.getNodes()[0],p.getNodes()[1]);
    sg.remove_link(-p.getNodes()[0],p.getNodes()[1]);
    sg.add_link(-src,new_node,old_link.dist);
    //TODO: migrate connection from the last middle node to the dest node into the new node
    old_link=sg.get_link(-p.getNodes()[p.getNodes().size()-2],p.getNodes()[p.getNodes().size()-1]);
    sg.remove_link(-p.getNodes()[p.getNodes().size()-2],p.getNodes()[p.getNodes().size()-1]);
    sg.add_link(-new_node,dest,old_link.dist);
    //TODO: delete nodes from the original path while there is any of them with a disconnecte end.
    bool mod=true;
    std::set<sgNodeID_t> mid_nodes;
    for (auto n:pmid.getNodes()) mid_nodes.insert(llabs(n));
    while (mod){
        mod=false;
        for (auto n:mid_nodes) {
            if (sg.get_fw_links(n).size()==0 or sg.get_bw_links(n).size()==0) {
                mod=true;
                sg.remove_node(n);
                edited_nodes.insert(llabs(n));
                mid_nodes.erase(n);
                break;
            }
        }
    }
}

void GraphEditor::join_path(SequenceGraphPath p, bool consume_nodes) {
    SequenceGraph& sg(ws.getGraph());
    std::set<sgNodeID_t> pnodes;
    for (auto n:p.getNodes()) {
        pnodes.insert( n );
        pnodes.insert( -n );
    }
    if (!p.is_canonical()) p.reverse();
    sgNodeID_t new_node=sg.add_node(Node(p.get_sequence()));
    //TODO:check, this may have a problem with a circle
    for (auto l:sg.get_bw_links(p.getNodes().front())) sg.add_link(new_node,l.dest,l.dist);

    for (auto l:sg.get_fw_links(p.getNodes().back())) sg.add_link(-new_node,l.dest,l.dist);

    //TODO: update read mappings
    if (consume_nodes) {
        for (auto n:p.getNodes()) {
            //check if the node has neighbours not included in the path.
            bool ext_neigh=false;
            if (n!=p.getNodes().back()) for (auto l:sg.get_fw_links(n)) if (pnodes.count(l.dest)==0) ext_neigh=true;
            if (n!=p.getNodes().front()) for (auto l:sg.get_bw_links(n)) if (pnodes.count(l.dest)==0) ext_neigh=true;
            if (ext_neigh) continue;
            sg.remove_node(n);
        }
    }
}