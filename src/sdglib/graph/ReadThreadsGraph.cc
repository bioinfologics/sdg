//
// Created by Bernardo Clavijo (EI) on 30/11/2020.
//

#include "ReadThreadsGraph.hpp"
#include <sdglib/views/NodeView.hpp>
#include <atomic>

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

NodeView ReadThreadsGraph::thread_start_nodeview(int64_t thread_id) {
    return get_nodeview((thread_id>0 ? thread_info[thread_id].start : thread_info[-thread_id].end));
}

NodeView ReadThreadsGraph::thread_end_nodeview(int64_t thread_id) {
    return thread_start_nodeview(-thread_id);
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

std::unordered_set<uint64_t> ReadThreadsGraph::node_threads(sgNodeID_t nid) {
    std::unordered_set<uint64_t> thread_ids;
    for (auto l: get_nodeview(nid).next()) thread_ids.insert(l.support().id);
    for (auto l: get_nodeview(nid).prev()) thread_ids.insert(l.support().id);
    return thread_ids;
}

std::vector<sgNodeID_t> ReadThreadsGraph::all_nids_fw_in_thread(sgNodeID_t nid, int64_t thread_id) {
    //TODO: this won't work if the thread has duplicated nodes.
    std::vector<sgNodeID_t> nodes;
    int ncf=0;
    int ncb=0;
    int np=-1;
    auto t=get_thread(thread_id);
    for (auto i=0;i<t.size();++i){
        if (t[i].node==nid) {
            ++ncf;
            np=i;
        }
        else if (t[i].node==nid) {
            ++ncb;
            np=i;
        }
    }
    if (ncf==1 and ncb==0){
        for (int i=np+1;i<t.size();++i) nodes.emplace_back(t[i].node);
    }
    if (ncf==0 and ncb==1){
        for (int i=np-1;i<-1;--i) nodes.emplace_back(-t[i].node);
    }
    return nodes;
}

ReadThreadsGraph ReadThreadsGraph::local_graph(sgNodeID_t nid, int64_t distance, uint16_t min_links) {
    ReadThreadsGraph lrtg(sdg,"local_rtg_"+std::to_string(nid)+"_"+std::to_string(distance)+"_"+std::to_string(min_links));
    nid=llabs(nid);
    for (auto tid:node_threads(nid)){
        auto thread=get_thread(tid);
        //find first start and last end of node
        int64_t first_start=10000000000,last_end=0;
        for (const auto &np:thread){
            if (llabs(np.node)==nid){
                if (np.start<first_start) first_start=np.start;
                if (np.end>last_end) last_end=np.end;
            }
        }
        //copy only nodes within the range into the new thread
        std::vector<NodePosition> new_thread;
        new_thread.reserve(thread.size());
        for (const auto &np:thread) {
            if (llabs(np.node)==nid or (np.start < first_start and np.end >= first_start - distance) or
                (np.end > last_end and np.start <= last_end + distance))
                new_thread.push_back(np);
        }
        //insert new thread into the local graph
        lrtg.add_thread(tid,new_thread);
    }
    //for each node in the new graph, if it doesnt have enough links, pop from all
    for (auto nv:lrtg.get_all_nodeviews(false,false)){
        if (lrtg.node_threads(nv.node_id()).size()<min_links) lrtg.pop_node_from_all(nv.node_id());
    }

    return lrtg;

}

std::vector<NodePosition> ReadThreadsGraph::get_thread(int64_t thread_id) {
    std::vector<NodePosition> thread;
    if (thread_info.count(llabs(thread_id))==0) return {};
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
            auto ln=next_in_thread(nid,llabs(thread_id),lc);
            nid=ln.node().node_id();
            p+=ln.distance();
            lc+=link_increment;
        }
        thread.emplace_back(nid,p,p+sdg.get_node_size(nid));
    } while (lc>=0 and lc<ti.link_count);
    return thread;
}

//TODO: this could be much faster by reversing the ends in thread info and updating the link index all along.
bool ReadThreadsGraph::flip_thread(int64_t thread_id) {
    thread_id=llabs(thread_id);
    if (thread_info.count(llabs(thread_id))==0) return false;
    auto rt=get_thread(-thread_id);
    remove_thread(thread_id);
    add_thread(thread_id,rt);
    return true;
}

std::unordered_map<sgNodeID_t, uint64_t> ReadThreadsGraph::node_thread_neighbours(sgNodeID_t nid) {
    std::unordered_map<sgNodeID_t, uint64_t> counts;
    for (auto tid:node_threads(nid)){
        for (auto &np:get_thread(tid)) counts[llabs(np.node)]+=1;
    }
    return counts;
}

int ReadThreadsGraph::clean_node(sgNodeID_t nid, int min_supported, int min_support) {
    auto n1tn=node_thread_neighbours(nid);
    std::vector<sgNodeID_t> to_pop;
    for (auto &tid:node_threads(nid)){
        int supported=0,unsupported=0;
        for (auto np:get_thread(tid)){
            if (n1tn[llabs(np.node)]>min_support) ++supported;
            else ++unsupported;
        }
        if (supported<min_supported) to_pop.emplace_back(tid);
    }

    for (auto tid:to_pop)
        pop_node(nid,tid);
    return to_pop.size();
}

std::vector<std::pair<uint64_t,sgNodeID_t>> ReadThreadsGraph::clean_all_nodes_popping_list( int min_supported,
                                                                                       int min_support) {
    std::vector<std::pair<uint64_t,sgNodeID_t>> popping_list;
    auto all_nvs=get_all_nodeviews(false,false);
    std::atomic<uint64_t> nc(0);
#pragma omp parallel
    {
        std::vector<std::pair<uint64_t,sgNodeID_t>> private_popping_list;
#pragma omp for
        for (auto i=0;i<all_nvs.size();++i){
            auto nid=all_nvs[i].node_id();
            auto n1tn=node_thread_neighbours(nid);
            for (auto &tid:node_threads(nid)){
                int supported=0,unsupported=0;
                for (auto np:get_thread(tid)){
                    if (n1tn[llabs(np.node)]>min_support) ++supported;
                    else ++unsupported;
                }
                if (supported<min_supported) private_popping_list.emplace_back(tid,nid);
            }
            if (++nc%10000) sdglib::OutputLog()<<nc<<" nodes analysed"<<std::endl;
        }
#pragma omp critical
        popping_list.insert(popping_list.end(),private_popping_list.cbegin(),private_popping_list.cend());

    };
    std::sort(popping_list.begin(),popping_list.end());
    return popping_list;
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
    auto thread_ids=node_threads(node_id);
    if (thread_ids.size()==0) return false;
    for(auto tid:thread_ids) pop_node(node_id,tid);
    return true;
}
