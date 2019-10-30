//
// Created by Bernardo Clavijo (EI) on 25/06/2018.
//

#include <functional>
#include "GraphEditor.hpp"

GraphEditorNodeExpansion::GraphEditorNodeExpansion(sgNodeID_t node,
                                                   std::vector<std::pair<sgNodeID_t, sgNodeID_t>> links_through) {

    if (node>0) {
        input_nodes.emplace_back(node);
    }
    else {
        input_nodes.emplace_back(-node);
    }
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

bool GraphEditor::detach_path(SequenceDistanceGraphPath p, bool consume_tips) {
    if (p.nodes.size()==1) return true;
    //TODO: check the path is valid?
    if (p.nodes.size()==2) {
        //std::cout<<"DETACHING path with only 2 nodes"<<std::endl;
        if (p.nodes.size()==2){
            auto fwls=ws.sdg.get_fw_links(p.nodes[0]);
            for (auto l:fwls) if (l.dest!=p.nodes[1]) ws.sdg.remove_link(l.source,l.dest);
            auto bwls=ws.sdg.get_bw_links(p.nodes[1]);
            for (auto l:bwls) if (l.dest!=-p.nodes[0]) ws.sdg.remove_link(l.source,l.dest);
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
    sgNodeID_t new_node=ws.sdg.add_node(Node(pmid.sequence()));
    //TODO: migrate connection from the src node to the first middle node into the new node
    auto old_link=ws.sdg.get_link(-p.nodes[0],p.nodes[1]);
    //ws.sdg.remove_link(-p.nodes[0],p.nodes[1]);
    auto old_fw_links=ws.sdg.get_fw_links(src);
    for (auto fwl:old_fw_links) ws.sdg.remove_link(fwl.source,fwl.dest);
    ws.sdg.add_link(-src,new_node,old_link.dist);
    //TODO: migrate connection from the last middle node to the dest node into the new node
    old_link=ws.sdg.get_link(-p.nodes[p.nodes.size()-2],p.nodes[p.nodes.size()-1]);
    //ws.sdg.remove_link(-p.nodes[p.nodes.size()-2],p.nodes[p.nodes.size()-1]);
    auto old_bw_links=ws.sdg.get_bw_links(dest);
    for (auto bwl:old_bw_links) ws.sdg.remove_link(bwl.source,bwl.dest);
    ws.sdg.add_link(-new_node,dest,old_link.dist);
    //TODO: delete nodes from the original path while there is any of them with a disconnecte end.
    bool mod=true;
    std::set<sgNodeID_t> mid_nodes;
    for (auto n:pmid.nodes) mid_nodes.insert(llabs(n));
    while (mod){
        mod=false;
        for (auto n:mid_nodes) {
            if (ws.sdg.get_fw_links(n).size()==0 or ws.sdg.get_bw_links(n).size()==0) {
                mod=true;
                ws.sdg.remove_node(n);
                edited_nodes.insert(llabs(n));
                mid_nodes.erase(n);
                break;
            }
        }
    }
    return true;
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
        //Create a unitig with the internal path sequence
        auto new_node=ws.sdg.add_node(SequenceDistanceGraphPath(ws.sdg,op.input_nodes).sequence());
        auto first_dist=ws.sdg.get_link(op.input_ends[0],op.input_nodes.front()).dist;
        auto last_dist=ws.sdg.get_link(op.input_ends[1],op.input_nodes.back()).dist;
        for (auto l:ws.sdg.get_fw_links(-op.input_ends[0])) ws.sdg.remove_link(l.source,l.dest);
        for (auto l:ws.sdg.get_bw_links(op.input_ends[1])) ws.sdg.remove_link(l.source,l.dest);
        ws.sdg.add_link(op.input_ends[0],new_node,first_dist);
        ws.sdg.add_link(-new_node,op.input_ends[1],last_dist);
    }
    node_expansion_queue.clear();
}

/*bool GraphEditor::detach_path(SequenceDistanceGraphPath p, bool consume_tips) {
    if (p.nodes.size()==1) return true;
    //TODO: check the path is valid?
    if (p.nodes.size()==2) {
        //std::cout<<"DETACHING path with only 2 nodes"<<std::endl;
        if (p.nodes.size()==2){
            auto fwls=ws.sdg.get_fw_links(p.nodes[0]);
            for (auto l:fwls) if (l.dest!=p.nodes[1]) ws.sdg.remove_link(l.source,l.dest);
            auto bwls=ws.sdg.get_bw_links(p.nodes[1]);
            for (auto l:bwls) if (l.dest!=-p.nodes[0]) ws.sdg.remove_link(l.source,l.dest);
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
    sgNodeID_t new_node=ws.sdg.add_node(Node(pmid.sequence()));
    //TODO: migrate connection from the src node to the first middle node into the new node
    auto old_link=ws.sdg.get_link(-p.nodes[0],p.nodes[1]);
    //ws.sdg.remove_link(-p.nodes[0],p.nodes[1]);
    auto old_fw_links=ws.sdg.get_fw_links(src);
    for (auto fwl:old_fw_links) ws.sdg.remove_link(fwl.source,fwl.dest);
    ws.sdg.add_link(-src,new_node,old_link.dist);
    //TODO: migrate connection from the last middle node to the dest node into the new node
    old_link=ws.sdg.get_link(-p.nodes[p.nodes.size()-2],p.nodes[p.nodes.size()-1]);
    //ws.sdg.remove_link(-p.nodes[p.nodes.size()-2],p.nodes[p.nodes.size()-1]);
    auto old_bw_links=ws.sdg.get_bw_links(dest);
    for (auto bwl:old_bw_links) ws.sdg.remove_link(bwl.source,bwl.dest);
    ws.sdg.add_link(-new_node,dest,old_link.dist);
    //TODO: delete nodes from the original path while there is any of them with a disconnecte end.
    bool mod=true;
    std::set<sgNodeID_t> mid_nodes;
    for (auto n:pmid.nodes) mid_nodes.insert(llabs(n));
    while (mod){
        mod=false;
        for (auto n:mid_nodes) {
            if (ws.sdg.get_fw_links(n).size()==0 or ws.sdg.get_bw_links(n).size()==0) {
                mod=true;
                ws.sdg.remove_node(n);
                edited_nodes.insert(llabs(n));
                mid_nodes.erase(n);
                break;
            }
        }
    }
    return true;
}*/

SequenceDistanceGraphPath GraphEditor::find_longest_path_from(sgNodeID_t node, std::string seq) {
    SequenceDistanceGraphPath p(ws.sdg);
    auto n1 = ws.sdg.nodes[llabs(node)];
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
        for (auto fwl:ws.sdg.get_fw_links(p.nodes.back())) {
            //Find the start position on the patch (use the connection distance)
            next_start=last_end+fwl.dist;
            //If the neighbour sequence and the patch sequence are the same till one of the two finishes: add neighbour as possible way forward.
            auto nn=ws.sdg.nodes[llabs(fwl.dest)];
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
            last_end=next_start+ws.sdg.nodes[llabs(next)].sequence.size();
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

void GraphEditor::join_path(SequenceDistanceGraphPath p, bool consume_nodes) {
    std::set<sgNodeID_t> pnodes;
    for (auto n:p.nodes) {
        pnodes.insert( n );
        pnodes.insert( -n );
    }
    if (!p.is_canonical()) p.reverse();
    sgNodeID_t new_node=ws.sdg.add_node(Node(p.sequence()));
    //TODO:check, this may have a problem with a circle
    for (auto l:ws.sdg.get_bw_links(p.nodes.front())) ws.sdg.add_link(new_node,l.dest,l.dist);

    for (auto l:ws.sdg.get_fw_links(p.nodes.back())) ws.sdg.add_link(-new_node,l.dest,l.dist);

    //TODO: update read mappings
    if (consume_nodes) {
        for (auto n:p.nodes) {
            //check if the node has neighbours not included in the path.
            bool ext_neigh=false;
            if (n!=p.nodes.back()) for (auto l:ws.sdg.get_fw_links(n)) if (pnodes.count(l.dest)==0) ext_neigh=true;
            if (n!=p.nodes.front()) for (auto l:ws.sdg.get_bw_links(n)) if (pnodes.count(l.dest)==0) ext_neigh=true;
            if (ext_neigh) continue;
            ws.sdg.remove_node(n);
        }
    }
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