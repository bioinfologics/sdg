//
// Created by Bernardo Clavijo (EI) on 30/11/2020.
//

#include "ReadThreadsGraph.hpp"
#include <sdglib/views/NodeView.hpp>

void ReadThreadsGraph::dump(std::string filename) {
    std::ofstream ofs(filename);
    DistanceGraph::write(ofs);
    sdglib::write_flat_unorderedmap(ofs,thread_info);
}

void ReadThreadsGraph::load(std::string filename) {
    std::ifstream ifs(filename);
    DistanceGraph::read(ifs);
    sdglib::read_flat_unorderedmap(ifs,thread_info);
}

bool ReadThreadsGraph::add_thread(int64_t thread_id, const std::vector<NodePosition> &node_positions, bool remove_duplicated, int min_thread_nodes) {
    std::unordered_set<sgNodeID_t> seen,duplicated;
    if (node_positions.size()<min_thread_nodes) return false;
    if (remove_duplicated) {
        seen.clear();
        duplicated.clear();
        for (const auto &p:node_positions){
            if (seen.count(llabs(p.node))) duplicated.insert(llabs(p.node));
            else seen.insert(llabs(p.node));
        }
    }
    uint16_t lidx=0;
    int16_t last_valid_i=-1;
    sgNodeID_t end1;
    for (auto i = 0; i < node_positions.size(); ++i) {
        if (remove_duplicated and  duplicated.count(llabs(node_positions[i].node))) continue;
        if (last_valid_i!=-1) {
            add_link(-node_positions[last_valid_i].node, node_positions[i].node,
                     node_positions[i].start - node_positions[last_valid_i].end,
                     {SupportType::LongRead, lidx++, thread_id});
        }
        else end1=node_positions[i].node;
        last_valid_i=i;
    }
    if (last_valid_i==-1) return false;
    sgNodeID_t end2=-node_positions[last_valid_i].node;
    thread_info[thread_id]={end1, end2, lidx};
    return true;
}

bool ReadThreadsGraph::remove_thread(int64_t thread_id) {
    thread_id=llabs(thread_id);
    if (thread_info.count(thread_id) == 0) return false;
    auto next_nid=thread_info[thread_id].start;
    int16_t next_link=0;
    bool link_deleted=false;
    do {
        link_deleted=false;
        for (auto nl:get_nodeview(next_nid).next()){
            if (nl.support().id==thread_id and nl.support().index==next_link) {
                remove_link(-next_nid,nl.node().node_id(),nl.distance(),nl.support());
                next_nid=nl.node().node_id();
                next_link+=1;
                link_deleted=true;
                break;
            }
        }

    } while (link_deleted);
    thread_info.erase(thread_id);
    return true;
}

NodeView ReadThreadsGraph::get_thread_start_nodeview(int64_t thread_id) {
    return get_nodeview((thread_id>0 ? thread_info[thread_id].start : thread_info[-thread_id].end));
}

NodeView ReadThreadsGraph::get_thread_end_nodeview(int64_t thread_id) {
    return get_thread_start_nodeview(-thread_id);
}

LinkView ReadThreadsGraph::next_in_thread(sgNodeID_t nid, int64_t thread_id, int64_t link_index) {
    bool found=false;
    sgNodeID_t nnid;
    Support s;
    int32_t d;
    for (auto l:links[llabs(nid)]) {
        if (l.source==-nid and l.support.id == thread_id and (link_index==-1 or l.support.index==link_index)) {
            if (found) throw std::runtime_error("There is more than one link out in thread");
            found=true;
            nnid=l.dest;
            d=l.dist;
            s=l.support;
        }
    }
    if (not found) throw std::runtime_error("There is no link out in thread");
    return LinkView(NodeView(this,nnid),d,s);
}

LinkView ReadThreadsGraph::prev_in_thread(sgNodeID_t nid, int64_t thread_id, int64_t link_index) {
    bool found=false;
    sgNodeID_t pnid;
    Support s;
    int32_t d;
    for (auto l:links[llabs(nid)]) {
        if (l.source==nid and l.support.id == thread_id and (link_index==-1 or l.support.index==link_index)) {
            if (found) throw std::runtime_error("There is more than one link out in thread");
            found=true;
            pnid=-l.dest;
            d=l.dist;
            s=l.support;
        }
    }
    if (not found) throw std::runtime_error("There is no link out in thread");
    return LinkView(NodeView(this,pnid),d,s);
}


std::vector<sgNodeID_t> ReadThreadsGraph::all_nids_fw_in_thread(sgNodeID_t nid, int64_t thread_id) {
    //TODO: this won't work if the thread has duplicated nodes.
    std::vector<sgNodeID_t> nodes;
    auto end1=-thread_info[llabs(thread_id)].start;
    auto end2=-thread_info[llabs(thread_id)].end;
    auto nnid=nid;
    while (nnid!=end1 and nnid!=end2) {
        nnid=next_in_thread(nid,thread_id).node().node_id();
        nodes.emplace_back(nnid);
    }
    return nodes;
}

std::vector<NodePosition> ReadThreadsGraph::get_thread(int64_t thread_id) {
    std::vector<NodePosition> thread;
    auto ti=thread_info[llabs(thread_id)];
    auto s=ti.start;
        int64_t lc=0;
    int link_increment=1;
    if (thread_id<0) {
        s=ti.end;
        link_increment=-1;
        lc=ti.link_count-1;
    }
    thread.reserve(ti.link_count);
    uint64_t p=0;
    sgNodeID_t nid=0;
    do {
        if (nid==0) nid=s;
        else {
            p+=sdg.get_node_size(nid);
            auto ln=next_in_thread(nid,thread_id,lc);
            nid=ln.node().node_id();
            p+=ln.distance();
            lc+=link_increment;
        }
        thread.emplace_back(nid,p,p+sdg.get_node_size(nid));
    } while (lc>=0 and lc<ti.link_count);
    return thread;
}

bool ReadThreadsGraph::pop_node(sgNodeID_t node_id, int64_t thread_id) {
    //Remove node, if node was start, update thread info to make the next node
    //check what happens when the node appears more than once in the thread
    //renumber links
    //update link count on threadinfo
    //or.... just get the thread out, remove the positions that have the node, and add the thread back.
    thread_id=llabs(thread_id);
    auto thread=get_thread(thread_id);
    std::vector<NodePosition> new_thread;
    new_thread.reserve(thread.size());
    for (auto np:thread) if (llabs(np.node)!=llabs(node_id)) new_thread.push_back(np);
    if (thread.size()==new_thread.size()) return false;
    remove_thread(thread_id);
    add_thread(thread_id,new_thread,false);
    return true;
}

bool ReadThreadsGraph::pop_node_from_all(sgNodeID_t node_id) {
    std::unordered_set<sgNodeID_t> thread_ids;
    for (auto l: get_nodeview(node_id).next()) thread_ids.insert(l.support().id);
    for (auto l: get_nodeview(node_id).prev()) thread_ids.insert(l.support().id);
    if (thread_ids.size()==0) return false;
    for(auto tid:thread_ids) pop_node(node_id,tid);
    return true;
}