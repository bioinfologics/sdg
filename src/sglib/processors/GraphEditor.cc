//
// Created by Bernardo Clavijo (EI) on 25/06/2018.
//

#include <functional>
#include "GraphEditor.hpp"

bool GraphEditor::detach_path(SequenceGraphPath p, bool consume_tips) {
    if (p.nodes.size()==1) return true;
    //TODO: check the path is valid?
    if (p.nodes.size()==2) {
        //std::cout<<"DETACHING path with only 2 nodes"<<std::endl;
        if (p.nodes.size()==2){
            auto fwls=ws.sg.get_fw_links(p.nodes[0]);
            for (auto l:fwls) if (l.dest!=p.nodes[1]) ws.sg.remove_link(l.source,l.dest);
            auto bwls=ws.sg.get_bw_links(p.nodes[1]);
            for (auto l:bwls) if (l.dest!=-p.nodes[0]) ws.sg.remove_link(l.source,l.dest);
        }
        return true;
    }
    //check if the path is already detached.
    if (p.is_unitig()) return true;
    //check if the nodes in the path have already been modified
    for (auto n:p.nodes) {
        if (edited_nodes.count(llabs(n))>0) {
            std::cout<<"NOT DETACHING: path has already edited-out nodes"<<std::endl;
            return false;
        }
    }
    //TODO: create a new node with the sequence of the "middle" nodes
    auto pmid=p;
    pmid.nodes.clear();
    for (auto i=1;i<p.nodes.size()-1;++i) pmid.nodes.push_back(p.nodes[i]);
    auto src=p.nodes.front();
    auto dest=p.nodes.back();
    if (!pmid.is_canonical()) {
        pmid.reverse();
        auto nsrc=-dest;
        dest=-src;
        src=nsrc;
    }
    sgNodeID_t new_node=ws.sg.add_node(Node(pmid.get_sequence()));
    //TODO: migrate connection from the src node to the first middle node into the new node
    auto old_link=ws.sg.get_link(-p.nodes[0],p.nodes[1]);
    ws.sg.remove_link(-p.nodes[0],p.nodes[1]);
    ws.sg.add_link(-src,new_node,old_link.dist);
    //TODO: migrate connection from the last middle node to the dest node into the new node
    old_link=ws.sg.get_link(-p.nodes[p.nodes.size()-2],p.nodes[p.nodes.size()-1]);
    ws.sg.remove_link(-p.nodes[p.nodes.size()-2],p.nodes[p.nodes.size()-1]);
    ws.sg.add_link(-new_node,dest,old_link.dist);
    //TODO: delete nodes from the original path while there is any of them with a disconnecte end.
    bool mod=true;
    std::set<sgNodeID_t> mid_nodes;
    for (auto n:pmid.nodes) mid_nodes.insert(llabs(n));
    while (mod){
        mod=false;
        for (auto n:mid_nodes) {
            if (ws.sg.get_fw_links(n).size()==0 or ws.sg.get_bw_links(n).size()==0) {
                mod=true;
                ws.sg.remove_node(n);
                edited_nodes.insert(llabs(n));
                mid_nodes.erase(n);
                break;
            }
        }
    }
}

int GraphEditor::patch_between(sgNodeID_t from, sgNodeID_t to, std::string patch) {
    auto n1 = ws.sg.nodes[llabs(from)];
    auto n2 = ws.sg.nodes[llabs(to)];
    const size_t ENDS_SIZE=1000;
    if (n1.sequence.size()>ENDS_SIZE) n1.sequence=n1.sequence.substr(n1.sequence.size()-ENDS_SIZE-1,ENDS_SIZE);
    if (n2.sequence.size()>ENDS_SIZE) n2.sequence.resize(ENDS_SIZE);
    if (from<0) n1.make_rc();
    if (to<0) n2.make_rc();
    auto p1=patch.find(n1.sequence);
    auto p2=patch.find(n2.sequence);
    if (p1>=patch.size() or p2>=patch.size() or p1>=p2) {
        return 1;
    }
    SequenceGraphPath sol(std::ref(ws.sg));
    auto paths=ws.sg.find_all_paths_between(from,to,patch.size(),30, false);
    //auto paths=ws.sg.find_all_paths_between(from,to,1000000,30);
    if (paths.empty()) return 2;
    for (auto p:paths){
        auto pnp=patch.find(p.get_sequence());
        if (pnp>p1 and pnp<p2) {
            if (sol.nodes.empty()) {
                sol.nodes.emplace_back(from);
                sol.nodes.insert(sol.nodes.end(),p.nodes.begin(),p.nodes.end());
                sol.nodes.emplace_back(to);
            }
            else {
                return 3;
            }
        }
    }
    if (sol.nodes.empty()) return 4;
    if (detach_path(sol,true)) return 0;
    else return 5;
}

void GraphEditor::join_path(SequenceGraphPath p, bool consume_nodes) {
    std::set<sgNodeID_t> pnodes;
    for (auto n:p.nodes) {
        pnodes.insert( n );
        pnodes.insert( -n );
    }
    if (!p.is_canonical()) p.reverse();
    sgNodeID_t new_node=ws.sg.add_node(Node(p.get_sequence()));
    //TODO:check, this may have a problem with a circle
    for (auto l:ws.sg.get_bw_links(p.nodes.front())) ws.sg.add_link(new_node,l.dest,l.dist);

    for (auto l:ws.sg.get_fw_links(p.nodes.back())) ws.sg.add_link(-new_node,l.dest,l.dist);

    //TODO: update read mappings
    if (consume_nodes) {
        for (auto n:p.nodes) {
            //check if the node has neighbours not included in the path.
            bool ext_neigh=false;
            if (n!=p.nodes.back()) for (auto l:ws.sg.get_fw_links(n)) if (pnodes.count(l.dest)==0) ext_neigh=true;
            if (n!=p.nodes.front()) for (auto l:ws.sg.get_bw_links(n)) if (pnodes.count(l.dest)==0) ext_neigh=true;
            if (ext_neigh) continue;
            ws.sg.remove_node(n);
        }
    }
}