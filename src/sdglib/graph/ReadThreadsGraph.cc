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
    int32_t last_valid_i=-1;
    sgNodeID_t end1;
    for (int32_t i = 0; i < node_positions.size(); ++i) {
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

NodeView ReadThreadsGraph::thread_start_nodeview(int64_t thread_id) const {
    return get_nodeview((thread_id>0 ? thread_info.at(thread_id).start : thread_info.at(-thread_id).end));
}

NodeView ReadThreadsGraph::thread_end_nodeview(int64_t thread_id) const {
    return thread_start_nodeview(-thread_id);
}

LinkView ReadThreadsGraph::next_in_thread(sgNodeID_t nid, int64_t thread_id, int64_t link_index) const {
    bool found=false;
    sgNodeID_t nnid;
    Support s;
    int32_t d;
    for (auto l:links[llabs(nid)]) {
        if (l.source==-nid and l.support.id == llabs(thread_id) and (link_index==-1 or l.support.index==link_index)) {
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

LinkView ReadThreadsGraph::prev_in_thread(sgNodeID_t nid, int64_t thread_id, int64_t link_index) const{
    bool found=false;
    sgNodeID_t pnid;
    Support s;
    int32_t d;
    for (auto l:links[llabs(nid)]) {
        if (l.source==nid and l.support.id == llabs(thread_id) and (link_index==-1 or l.support.index==link_index)) {
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

std::unordered_set<int64_t> ReadThreadsGraph::node_threads(sgNodeID_t nid, bool oriented) const {
    std::unordered_set<int64_t> thread_ids;
    for (const auto &lc:{get_nodeview(nid).next(),get_nodeview(nid).prev()}) {
        for (const auto &l:lc) {
            const auto &t = l.support().id;
            if (oriented) thread_ids.insert(thread_fw_in_node(t,nid) ? t:-t);
            else thread_ids.insert(t);
        }
    }
    return thread_ids;
}

std::unordered_set<std::pair<int64_t, int32_t>> ReadThreadsGraph::node_threadpositions(sgNodeID_t nid) const {
    std::unordered_map<int64_t,std::pair<int32_t,int32_t>> thread_links; //Populate a map with all threads and their positions left and right of this node
    thread_links.reserve(links[llabs(nid)].size());
    for (const auto &l:links[llabs(nid)]) {
        if (l.source==nid) thread_links[l.support.id].first=l.support.index+1;//XXX: This is to avoid problems because links are counted from 0
        else thread_links[l.support.id].second=l.support.index+1;
    }
    std::unordered_set<std::pair<int64_t, int32_t>> thread_positions;
    for (const auto &tl:thread_links){
        if (tl.second.first==0) {
            if (tl.second.second==1) thread_positions.emplace(tl.first,1);
            else thread_positions.emplace(-tl.first,1);
        }
        else if (tl.second.second==0) {
            if (tl.second.first==1) thread_positions.emplace(-tl.first,thread_info.at(llabs(tl.first)).link_count+1);
            else thread_positions.emplace(tl.first,thread_info.at(llabs(tl.first)).link_count+1);
        }
        else if (tl.second.first<tl.second.second) thread_positions.emplace(tl.first,tl.second.second);
        else thread_positions.emplace(-tl.first,thread_info.at(llabs(tl.first)).link_count+1-tl.second.second);
    }
    return thread_positions;
}

std::vector<sgNodeID_t> ReadThreadsGraph::all_nids_fw_in_thread(sgNodeID_t nid, int64_t thread_id) {
    //TODO: this won't work if the thread has duplicated nodes.
    std::vector<sgNodeID_t> nodes;
    int ncf=0;
    int ncb=0;
    int np=-1;
    auto t=get_thread(thread_id);
    for (auto i=0;i<t.size();++i){
        if (t[i].node==-nid) {
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

std::vector<NodePosition> ReadThreadsGraph::get_thread(int64_t thread_id) const{
    std::vector<NodePosition> thread;
    if (thread_info.count(llabs(thread_id))==0) return {};
    auto ti=thread_info.at(llabs(thread_id));
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

std::unordered_map<sgNodeID_t, uint64_t> ReadThreadsGraph::node_thread_neighbours(sgNodeID_t nid, bool oriented) {
    std::unordered_map<sgNodeID_t, uint64_t> counts;
    for (auto tid:node_threads(nid)){
        if (not oriented) {
            for (auto &np:get_thread(tid)) ++counts[llabs(np.node)];
        }
        else {
            auto t=get_thread(tid);
            bool fw=true;
            for (auto &np:t) {
                if (np.node==-nid){
                    fw=false;
                    break;
                }
                else if (np.node==nid) break;
            }
            for (auto &np:t) ++counts[(fw ? np.node:-np.node)];
        }
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

std::vector<std::pair<uint64_t, sgNodeID_t>> ReadThreadsGraph::clean_repeat_nodes_popping_list(int max_threads) {
    std::vector<std::pair<uint64_t,sgNodeID_t>> popping_list;
    auto all_nvs=get_all_nodeviews(false,false);
    std::atomic<uint64_t> nc(0);
#pragma omp parallel
    {
        std::vector<std::pair<uint64_t,sgNodeID_t>> private_popping_list;
#pragma omp for
        for (auto i=0;i<all_nvs.size();++i){
            auto nid=all_nvs[i].node_id();
            auto nts=node_threads(nid);
            if (nts.size()>max_threads){
                for (auto tid:nts)
                    private_popping_list.emplace_back(tid,nid);
            }
            if (++nc%10000==0) sdglib::OutputLog()<<nc<<" nodes analysed"<<std::endl;
        }
#pragma omp critical
        popping_list.insert(popping_list.end(),private_popping_list.cbegin(),private_popping_list.cend());

    };
    std::sort(popping_list.begin(),popping_list.end());
    return popping_list;
}

std::vector<std::pair<uint64_t, sgNodeID_t>> ReadThreadsGraph::clean_lowmapping_nodes_popping_list(int min_threads) {
    std::vector<std::pair<uint64_t,sgNodeID_t>> popping_list;
    auto all_nvs=get_all_nodeviews(false,false);
    std::atomic<uint64_t> nc(0);
#pragma omp parallel
    {
        std::vector<std::pair<uint64_t,sgNodeID_t>> private_popping_list;
#pragma omp for
        for (auto i=0;i<all_nvs.size();++i){
            auto nid=all_nvs[i].node_id();
            auto nts=node_threads(nid);
            if (nts.size()<min_threads){
                for (auto tid:nts)
                    private_popping_list.emplace_back(tid,nid);
            }
            if (++nc%10000==0) sdglib::OutputLog()<<nc<<" nodes analysed"<<std::endl;
        }
#pragma omp critical
        popping_list.insert(popping_list.end(),private_popping_list.cbegin(),private_popping_list.cend());

    };
    std::sort(popping_list.begin(),popping_list.end());
    return popping_list;
}

std::vector<std::pair<int64_t,sgNodeID_t>> ReadThreadsGraph::clean_all_nodes_popping_list( int min_supported,
                                                                                       int min_support) {
    std::vector<std::pair<int64_t,sgNodeID_t>> popping_list;
    auto all_nvs=get_all_nodeviews(false,false);
    std::atomic<uint64_t> nc(0);
#pragma omp parallel
    {
        std::vector<std::pair<int64_t,sgNodeID_t>> private_popping_list;
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
                if (supported<=min_supported) private_popping_list.emplace_back(tid,nid);
            }
            if (++nc%10000==0) sdglib::OutputLog()<<nc<<" nodes analysed"<<std::endl;
        }
#pragma omp critical
        popping_list.insert(popping_list.end(),private_popping_list.cbegin(),private_popping_list.cend());

    };
    std::sort(popping_list.begin(),popping_list.end());
    return popping_list;
}

std::unordered_map<uint64_t,std::set<sgNodeID_t>> ReadThreadsGraph::thread_nodesets(){
    std::unordered_map<uint64_t,std::set<sgNodeID_t>> tns;
    for (auto &nv:get_all_nodeviews(false,false)) {
        auto nid = nv.node_id();
        for (auto &tid:node_threads(nid)) {
            tns[tid].insert(nid);
        }
    }
    return tns;
}

size_t nodeset_intersection_size(const std::set<sgNodeID_t>& v1, const std::set<sgNodeID_t>& v2)
{
    size_t s=0;
    for (auto p1=v1.cbegin(),p2=v2.cbegin();p1!=v1.cend() and p2!=v2.cend();){
        if (*p1==*p2) {
            ++s;
            ++p1;
            ++p2;
        }
        else if (*p1<*p2) ++p1;
        else ++p2;
    }
    return s;
}

std::vector<std::pair<int64_t, sgNodeID_t>> ReadThreadsGraph::clean_all_nodes_by_thread_clustering_popping_list(
        int min_shared, float max_second_perc) {
    std::vector<std::pair<int64_t,sgNodeID_t>> popping_list;
    auto all_nvs=get_all_nodeviews(false,false);
    std::atomic<uint64_t> nc(0);
    std::atomic<uint64_t> second_too_close(0);
    std::atomic<uint64_t> second_second_far_enough(0);
    auto tns=thread_nodesets();
#pragma omp parallel
    {
        std::vector<std::pair<int64_t,sgNodeID_t>> private_popping_list;
        std::unordered_map<uint64_t,std::vector<uint64_t>> tconns;
#pragma omp for
        for (auto i=0;i<all_nvs.size();++i){
            auto nid=all_nvs[i].node_id();
            auto ntids=node_threads(nid);

            //first analyse connections between threads (all-vs-all, but restricted to the node, caching would work too)
            tconns.clear();
            for (auto &tid1:ntids)
                for (auto &tid2:ntids)
                    if (nodeset_intersection_size(tns[tid1],tns[tid2])>=min_shared) tconns[tid1].emplace_back(tid2);

            //now cluster threads in the more rudimentary way possible
            std::vector<std::set<int64_t>> tclusters;
            std::unordered_set<uint64_t> used_tids;
            std::vector<uint64_t>to_explore_tids;
            for (auto &tid1:ntids){
                if (used_tids.count(tid1)) continue;
                tclusters.push_back({tid1});
                used_tids.insert(tid1);
                auto to_explore_tids=tconns[tid1];
                while (not to_explore_tids.empty()){
                    std::vector<uint64_t> next_to_explore;
                    for (auto &tid2:to_explore_tids){
                        if (used_tids.count(tid2)) continue;
                        tclusters.back().insert(tid2);
                        used_tids.insert(tid2);
                        for (auto &tid3:tconns[tid2]){
                            if (used_tids.count(tid3)) continue;
                            next_to_explore.push_back(tid3);
                        }
                    }
                    to_explore_tids=next_to_explore;
                }
            }
            //now analyse this clustering
            if (tclusters.size()==1) continue; //all is as should be
            //reverse ordering by size
            std::sort(tclusters.begin(),tclusters.end(),[](const std::set<int64_t> &a, const std::set<int64_t> &b){return a.size() > b.size();});
            //now take those decisions
            //std::cout << "Cluster size comparison: " << tclusters[0].size() << " * " << max_second_perc << " < " << tclusters[1].size() << std::endl;
            if (tclusters[0].size()*max_second_perc<tclusters[1].size()){
                //second cluster is too big, pop from all!
                second_too_close++;
                for (auto &tid:ntids){
                    private_popping_list.push_back({tid,nid});
                }
            }
            else {
                //pop from everywhere but the first cluster
                second_second_far_enough++;
                for (auto &tid:ntids){
                    if (tclusters[0].count(tid)) continue;
                    private_popping_list.push_back({tid,nid});
                }
            }

            if (++nc%10000==0) sdglib::OutputLog()<<nc<<" nodes analysed"<<std::endl;
        }
#pragma omp critical
        popping_list.insert(popping_list.end(),private_popping_list.cbegin(),private_popping_list.cend());

    };
    std::sort(popping_list.begin(),popping_list.end());
    std::cout << "Second far enough: " << second_second_far_enough << ", second too close: " << second_too_close <<std::endl;
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

bool ReadThreadsGraph::pop_nodes(std::vector<sgNodeID_t> node_ids, int64_t thread_id) {
    //Remove node, if node was start, update thread info to make the next node
    //check what happens when the node appears more than once in the thread
    //renumber links
    //update link count on threadinfo
    //or.... just get the thread out, remove the positions that have the node, and add the thread back.
    std::unordered_set<sgNodeID_t> nids;
    for (auto &nid:node_ids) nids.insert(llabs(nid));
    thread_id=llabs(thread_id);
    auto thread=get_thread(thread_id);
    std::vector<NodePosition> new_thread;
    new_thread.reserve(thread.size());
    for (auto np:thread) if (nids.count(llabs(np.node))==0) new_thread.push_back(np);
    if (thread.size()==new_thread.size()) return false;
    remove_thread(thread_id);
    add_thread(thread_id,new_thread,false);
    return true;
}

void ReadThreadsGraph::apply_popping_list(const std::vector<std::pair<uint64_t, sgNodeID_t>> &popping_list) {
    //if list is not sorted, this will be as slow as popping each node, or more!
    std::vector<sgNodeID_t> node_ids;
    if (popping_list.empty()) return;
    int64_t last_tid=popping_list[0].first;
    for (auto &tn:popping_list){
        if (tn.first!=last_tid){
            pop_nodes(node_ids,last_tid);
            node_ids.clear();
            last_tid=tn.first;
        }
        node_ids.emplace_back(tn.second);

    }
    pop_nodes(node_ids,last_tid);
}

bool ReadThreadsGraph::pop_node_from_all(sgNodeID_t node_id) {
    auto thread_ids=node_threads(node_id);
    if (thread_ids.size()==0) return false;
    for(auto tid:thread_ids) pop_node(node_id,tid);
    return true;
}

bool ReadThreadsGraph::thread_fw_in_node(int64_t tid, sgNodeID_t nid) const {
    if (llabs(nid)>links.size()) throw std::runtime_error("thread "+std::to_string(tid)+" not in node "+std::to_string(nid));
    if (thread_info.at(llabs(tid)).link_count==1) {//for single link threads, we need to evaluate the orientation from the ends
        if (tid>0) return nid==thread_info.at(llabs(tid)).start or nid==-thread_info.at(llabs(tid)).end;
        return nid==-thread_info.at(llabs(tid)).start or nid==thread_info.at(llabs(tid)).end;
    }
    int prev_idx=-1;
    int next_idx=-1;
    auto atid=llabs(tid);
    for (auto &l:links[llabs(nid)]) {
        if (l.support.id==atid) {
            if (l.source==nid) prev_idx=l.support.index;
            else next_idx=l.support.index;
        }
    }
    if (prev_idx==-1 and next_idx==-1) throw std::runtime_error("thread "+std::to_string(tid)+" not in node "+std::to_string(nid));
    if (prev_idx==-1 and next_idx==0) return tid>0;
    if (prev_idx==-1 and next_idx!=0) return tid<0;
    if (prev_idx==0 and next_idx==-1) return tid<0;
    if (prev_idx!=0 and next_idx==-1) return tid>0;
    if (prev_idx<next_idx) return tid>0;
    return tid<0;

}

int64_t ReadThreadsGraph::nodes_before_in_thread(int64_t tid, sgNodeID_t nid) const {
    if (llabs(nid)>links.size()) throw std::runtime_error("thread "+std::to_string(tid)+" not in node "+std::to_string(nid));
    int prev_idx=-1;
    int next_idx=-1;
    auto atid=llabs(tid);
    for (auto &l:links[llabs(nid)]) {
        if (l.support.id==atid) {
            if (l.source==nid) prev_idx=l.support.index;
            else next_idx=l.support.index;
        }
    }
    auto thread_links=thread_info.at(llabs(tid)).link_count;
    if (prev_idx==-1 and next_idx!=-1) return 0;
    if (next_idx==-1 and prev_idx!=-1) return thread_links;
    if (prev_idx<next_idx and tid>0) return prev_idx+1;
    if (prev_idx>next_idx and tid<0) return thread_links-prev_idx;
    throw std::runtime_error("node is not in thread!");
}

int64_t ReadThreadsGraph::nodes_after_in_thread(int64_t tid, sgNodeID_t nid) const {
    if (llabs(nid)>links.size()) throw std::runtime_error("thread "+std::to_string(tid)+" not in node "+std::to_string(nid));
    int prev_idx=-1;
    int next_idx=-1;
    auto atid=llabs(tid);
    for (auto &l:links[llabs(nid)]) {
        if (l.support.id==atid) {
            if (l.source==nid) prev_idx=l.support.index;
            else next_idx=l.support.index;
        }
    }
    auto thread_links=thread_info.at(llabs(tid)).link_count;
    if (prev_idx!=-1 and next_idx==-1) return 0;
    if (next_idx!=-1 and prev_idx==-1) return thread_links;
    if (prev_idx<next_idx and tid>0) return thread_links-next_idx;
    if (prev_idx>next_idx and tid<0) return prev_idx;
    throw std::runtime_error("node is not in thread!");
}

//TODO:: split into make_thread_nodepositions / pick_next (which automatically calls clean if needed) / -> make it a small class??
std::map<uint64_t, std::vector<std::pair<int64_t, sgNodeID_t>>> ReadThreadsGraph::make_thread_nodepositions(const std::set<sgNodeID_t> & nodes) const{
    std::map<uint64_t, std::vector<std::pair<int64_t, sgNodeID_t>>> thread_node_positions;
    //1- for every thread id find the lowest and the highest link number asociated with any node.
    std::map<uint64_t ,std::pair<int64_t,int64_t>> thread_limits;
    for (auto nid:nodes){
        auto nv=get_nodeview(nid);
        // Of all nodes, links fw and bw, reads ids are the map kays, map values are the first and last index
        // at the end map contains the first and lat index in the read map to a node (not specified wich node)
        for (auto &lv: {nv.next(),nv.prev()}) {
            for (auto &l: lv) {
                if (thread_limits.count(l.support().id) == 0)
                    thread_limits[l.support().id] = {l.support().index, l.support().index};
                else if (thread_limits[l.support().id].first > l.support().index)
                    thread_limits[l.support().id].first = l.support().index;
                else if (thread_limits[l.support().id].second < l.support().index)
                    thread_limits[l.support().id].second = l.support().index;
            }
        }
    }
    // Delete limit pairs that are consecutive or the same number
    for (auto it = thread_limits.cbegin(); it != thread_limits.cend(); ){
        if (it->second.second<=it->second.first+1) it=thread_limits.erase(it++);
        else ++it;
    }
    //TODO: count all fw/bw in node for each thread and discard those where nodes appear in more than one direciton.
//    for (auto &tl:thread_limits){
//        std::cout<<"Thread "<<tl.first<<" -> "<<tl.second.first<<" : "<<tl.second.second<<std::endl;
//    }
    //Find the first nv with a next in each thread -> prev  is limit/outside/inexistant
    for (auto nid:nodes) {
        auto nv=get_nodeview(nid);
        for (auto ln:nv.next()){
            auto tl=thread_limits.find(ln.support().id);
            if (tl==thread_limits.end()) continue;
            auto tp=thread_node_positions.find(ln.support().id);
            if (thread_fw_in_node(ln.support().id,nid)) {
                if ( (tp==thread_node_positions.end() and ln.support().index==tl->second.first+1) or ln.support().index==tl->second.first)
                    thread_node_positions[ln.support().id]={{0,nid}};
            }
            else {
                if ( (tp==thread_node_positions.end() and ln.support().index==tl->second.second-1) or ln.support().index==tl->second.second)
                    thread_node_positions[ln.support().id]={{0,nid}};
            }
        }

    }

    for (auto &tnp:thread_node_positions){
        auto &tl=thread_limits[tnp.first];
        int last_end_p=sdg.get_node_size(tnp.second[0].second);
        auto ln=next_in_thread(tnp.second.back().second,tnp.first);
        while (ln.support().index>=tl.first and ln.support().index<=tl.second){
            if (nodes.count(ln.node().node_id())) {
                tnp.second.emplace_back(last_end_p+ln.distance(),ln.node().node_id());
            }

            last_end_p+=ln.distance()+ln.node().size();

            try {
                ln=next_in_thread(ln.node().node_id(),tnp.first);
            }
            catch (const std::exception&) {
                break;
            }
        }
    }
    return thread_node_positions;

}