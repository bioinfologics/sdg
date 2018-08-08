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
    //ws.sg.remove_link(-p.nodes[0],p.nodes[1]);
    auto old_fw_links=ws.sg.get_fw_links(src);
    for (auto fwl:old_fw_links) ws.sg.remove_link(fwl.source,fwl.dest);
    ws.sg.add_link(-src,new_node,old_link.dist);
    //TODO: migrate connection from the last middle node to the dest node into the new node
    old_link=ws.sg.get_link(-p.nodes[p.nodes.size()-2],p.nodes[p.nodes.size()-1]);
    //ws.sg.remove_link(-p.nodes[p.nodes.size()-2],p.nodes[p.nodes.size()-1]);
    auto old_bw_links=ws.sg.get_bw_links(dest);
    for (auto bwl:old_bw_links) ws.sg.remove_link(bwl.source,bwl.dest);
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
    return true;
}

SequenceGraphPath GraphEditor::find_longest_path_from(sgNodeID_t node, std::string seq) {
    SequenceGraphPath p(ws.sg);
    auto n1 = ws.sg.nodes[llabs(node)];
    const size_t ENDS_SIZE=200;
    if (node<0) n1.make_rc();
    if (n1.sequence.size()>ENDS_SIZE) n1.sequence=n1.sequence.substr(n1.sequence.size()-ENDS_SIZE,ENDS_SIZE);

    //Find the end point of the current node in the patch sequence
    auto p1=seq.find(n1.sequence);
    if (p1>=seq.size()) return p;
    //std::cout<<"Found "<<n1.sequence.size()<<" bp of starting node at position "<<p1<<" in the patch"<<std::endl;
    auto last_end=p1+n1.sequence.size();
    //std::cout<<"last_end is now "<<last_end<<std::endl;
    //std::cout<<"n1 seq:    "<<n1.sequence<<std::endl;
    //std::cout<<"patch seq: "<<seq<<std::endl;
    p.nodes.emplace_back(node);

    while (last_end<seq.size()) {
        sgNodeID_t next=0;
        uint64_t next_start;
        //For every fw connection:
        for (auto fwl:ws.sg.get_fw_links(p.nodes.back())) {
            //Find the start position on the patch (use the connection distance)
            next_start=last_end+fwl.dist;
            //If the neighbour sequence and the patch sequence are the same till one of the two finishes: add neighbour as possible way forward.
            auto nn=ws.sg.nodes[llabs(fwl.dest)];
            if (fwl.dest<0) nn.make_rc();
            auto s=seq.c_str()+next_start;
            auto ns=nn.sequence.c_str();
            //std::cout<<"Evaluating possible next node: "<<fwl.dest<<" distance is "<< fwl.dist<<" which makes next_start="<<next_start<<std::endl;
            //std::cout<<"Node seq:  "<<ns<<std::endl;
            //std::cout<<"Patch seq: "<<s<<std::endl;
            while (*s!='\0' and *ns!='\0' and *s==*ns) {++s;++ns;};
            if (*s=='\0' or *ns=='\0') {
                if (next==0) next=fwl.dest;
                else {
                    next=0;
                    break;
                }
            }
        }
        if (next==0) break;
        else {
            last_end=next_start+ws.sg.nodes[llabs(next)].sequence.size();
            p.nodes.emplace_back(next);
        }
    }
    //return path
    return p;
}

int GraphEditor::patch_between(sgNodeID_t from, sgNodeID_t to, std::string patch) {
    auto p=find_longest_path_from(from,patch);
    //std::cout<<"trying to patch between "<<from<<" and "<< to<<", longest patch:";
    //for (auto n:p.nodes) std::cout<<" "<<n;
    //std::cout<<std::endl;
    if (p.nodes.empty()) return 1;
    if (p.nodes.back()!=to) return 2;
    if (detach_path(p,true)) return 0;
    else return 5;
};

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