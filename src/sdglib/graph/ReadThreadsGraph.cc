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

std::map<sgNodeID_t,std::pair<uint64_t,uint64_t>> ReadThreadsGraph::make_node_first_later(const std::map<uint64_t,std::vector<std::pair<int64_t,sgNodeID_t>>> &thread_node_positions, const std::map<uint64_t,int64_t> &thread_nextpos) const{
    std::map<sgNodeID_t,std::pair<uint64_t,uint64_t>> node_first_later;
    for (const auto &tnp:thread_node_positions){
        int start=0;
        if (thread_nextpos.count(tnp.first)) start=thread_nextpos.at(tnp.first);
        if (start==-1 or start>=tnp.second.size()) continue;

        ++node_first_later[tnp.second[start].second].first;
        for (auto i=start+1;i<tnp.second.size();++i) ++node_first_later[tnp.second[i].second].second;
    }
    return node_first_later;
}

bool ReadThreadsGraph::clean_thread_nodepositions(std::map<uint64_t,std::vector<std::pair<int64_t,sgNodeID_t>>> &thread_node_positions,std::set<sgNodeID_t> nodes_to_review) const {
    bool cleaned = false;
    while (true) {
//        std::cout << "cleaning positions in " << thread_node_positions.size() << " threads for "
//                  << nodes_to_review.size() << " nodes" << std::endl;
        std::map<sgNodeID_t, std::pair<std::map<sgNodeID_t, int>, std::map<sgNodeID_t, int>>> node_prev_next;

        //initialise total counts aggregating over each thread
        for (auto &tnp:thread_node_positions) {
//            std::cout << "Thread " << tnp.first << ": ";
            for (auto &pn:tnp.second) {
//                std::cout << pn.second << ", ";
                if (nodes_to_review.count(pn.second)) {
                    bool before = true;
                    for (auto &n:tnp.second) {
                        if (before) {
                            if (n != pn) ++node_prev_next[pn.second].first[n.second];
                            else before = false;
                        } else ++node_prev_next[pn.second].second[n.second];
                    }
                }
            }
//            std::cout << std::endl;
        }

//    for (auto &npn:node_prev_next){
//        for (auto &np:npn.second.first) if (npn.second.second.count(np.first)) std::cout<<"Potential conflict for node "<<npn.first<<" "<<np.second<<" : " << npn.second.second[np.first]<<std::endl;
//    }

        //initialise total counts aggregating over each thread
        int64_t thread = 0;
        sgNodeID_t node = 0;
        int max_conflicts = 0;
        for (auto &tnp:thread_node_positions) {

            for (auto &pn:tnp.second) {
                int conflicts = 0;
                if (nodes_to_review.count(pn.second)) {
                    bool before = true;
                    for (auto &n:tnp.second) {
                        if (before) {
                            if (n == pn) before = false;
                            else if (node_prev_next[pn.second].second[n.second]) conflicts+=node_prev_next[pn.second].second[n.second];
                        } else if (node_prev_next[pn.second].first[n.second]) conflicts+=node_prev_next[pn.second].first[n.second];
                    }
                }
                if (conflicts > max_conflicts) {
//                    std::cout << "Node " << pn.second << " has " << conflicts << " conflicts vs. >1 thread on thread "
//                              << tnp.first << std::endl;
                    thread = tnp.first;
                    node = pn.second;
                    max_conflicts = conflicts;
                }
            }
        }
        if (max_conflicts > 0) {
            //remove node from thread
//            std::cout<<"removing node "<<node<<" from thread "<<thread<<" with "<<max_conflicts<<" conflicts"<<std::endl;
            thread_node_positions[thread].erase(std::remove_if(thread_node_positions[thread].begin(),thread_node_positions[thread].end(),[&node](auto &e){return e.second==node;}));
            cleaned=true;
        }
        else break;
    }
    std::vector<int64_t> to_delete;
    for (auto &tnp:thread_node_positions) {
        if (tnp.second.size()<2) to_delete.push_back(tnp.first);
    }
    for (auto &tid:to_delete) thread_node_positions.erase(tid);
    return cleaned;
}

std::vector<sgNodeID_t> ReadThreadsGraph::order_nodes(const std::vector<sgNodeID_t> nodes, bool write_detailed_log) const {
    //rustic attempt at ordering the nodes for now, but using every read between them.

    std::set<sgNodeID_t> nodeset(nodes.begin(),nodes.end());

    auto thread_node_positions=make_thread_nodepositions(nodeset);

    auto node_first_later=make_node_first_later(thread_node_positions);
    clean_thread_nodepositions(thread_node_positions,nodeset);
    node_first_later=make_node_first_later(thread_node_positions);
    std::map<uint64_t,int64_t> thread_nextpos;
    for (auto &tnp:thread_node_positions) thread_nextpos[tnp.first]=0;
    std::ofstream log;
    if (write_detailed_log) {
        log.open("order_nodes_log.txt",std::ios_base::out|std::ios_base::app);
        log<<std::endl<<"---"<<std::endl<<"Order nodes: [";
        for (auto &n:nodes) log<<n<<", ";
        log<<"]"<<std::endl;
    }


    std::vector<sgNodeID_t> sorted_nodes;
    std::set<sgNodeID_t> candidates;

    //Just iterate finding the next node to put in:
    while (true){
        sgNodeID_t nid=0; // basically, once this is found, we're done for the iteration

        //check if all threads are finished already
        bool finished=true;
        for (auto &tnp:thread_nextpos) if (tnp.second!=-1) finished=false;
        if (finished){
            if (write_detailed_log==true) log<<"All threads finished, order done."<<std::endl;
            break;
        }

        //first look at nodes that only appear at the front of their threads.
        candidates.clear();
        for (auto &nfl:node_first_later) if (nfl.second.first and not nfl.second.second) candidates.insert(nfl.first);

        // TRIVIAL CASE: just one node is at the front only.
        if (candidates.size()==1){
            nid=*candidates.begin();
            log<<"T "<<nid<<std::endl;
        }

        // DISTANCE ANALYSIS REQUIRED: more than one candidate at the front only
        else if (candidates.size()>1){ //Multiple candidates untie by distance.
            if (write_detailed_log) {
                log << "M [";
                for (auto &n:candidates) log << n << ", ";
                log<<"] ";
            }
            //find a prevs that are connected to all candidates, or nexts same. pick first by distance. all picks must agree
            std::map<std::pair<sgNodeID_t ,sgNodeID_t>,std::vector<int64_t> > dists; //<node,prev>->count;
            for (auto &tn:thread_nextpos){
                auto tnnid=thread_node_positions[tn.first][tn.second].second;
                if (tn.second!=-1 and candidates.count(tnnid)){
                    auto tnpos=thread_node_positions[tn.first][tn.second].first;
                    for (auto &onp:thread_node_positions[tn.first])
                        if (onp.second!=tnnid) dists[{tnnid,onp.second}].emplace_back(onp.first-tnpos);
                }
            }
            std::map<sgNodeID_t, uint64_t> dcounts;
            for (auto &d:dists) ++dcounts[d.first.second];

            //this is the all picks must agree part
            std::map<sgNodeID_t ,int> winner_scores;
            sgNodeID_t winner=0;
            for (auto oc:dcounts){
                auto local_winner=0;
                if (oc.second!=candidates.size()) continue; //not in all candidates
                int64_t largest_dist=INT64_MIN;
                int64_t min_support=INT64_MAX;
                for (auto c:candidates){
                    //poor man's median
                    auto &dv=dists[{c,oc.first}];
                    std::sort(dv.begin(),dv.end());
                    auto d=dv[dv.size()/2];
                    min_support=std::min(min_support,(int64_t) dv.size());
                    if (d>largest_dist){
                        largest_dist=d;
                        local_winner=c;
                    }
                }
                /*if (winner==0) winner=local_winner;
                if (winner!=local_winner) {
                    std::cout<<"aborting order, untie has different winners!"<<std::endl;
                    winner=0;
                    break;
                }*/
                winner_scores[local_winner]+=min_support;

            }
            for (auto &ws:winner_scores) if (winner==0 or ws.second>winner_scores[winner]) winner=ws.first;
            for (auto &ws:winner_scores) {
                if (ws.second == winner_scores[winner] and ws.first != winner) {
                    if (write_detailed_log) log<<"aborting order, untie has two different winners!"<<std::endl;
                    winner = 0;
                    break;
                }
            }
            if (winner==0) {
                if (write_detailed_log) log<<"aborting order, can't find a winner!"<<std::endl;
                break;
            }
            log<<winner<<std::endl;
            nid=winner;
        }

        // THREAD DISORDER: no candidate appears front-only, threads must be filtered to be order-coherent.
        else if (candidates.size()==0){
            if (write_detailed_log) log<<"C"<<std::endl;


            //TODO: just call clean_threadnodepositions from here.
//            std::cout<<"Nodes to order: ";
//            for (auto &n:nodes) std::cout<<n<<", ";
//            std::cout<<std::endl;
//            std::cout<<"Order so far: ";
//            for (auto &n:sorted_nodes) std::cout<<n<<", ";
//            std::cout<<std::endl<<std::endl;
//            for (auto &nfl:node_first_later) std::cout<<"Node "<<nfl.first<<": f="<<nfl.second.first<<" l="<<nfl.second.second<<std::endl;
//            std::cout<<std::endl;
//            for (auto &tnp:thread_node_positions) {
//                std::cout<<"Thread "<<tnp.first<<": ";
//                for (auto i=thread_nextpos[tnp.first];i<tnp.second.size();++i)std::cout<<tnp.second[i].second<<", ";
//                std::cout<<std::endl;
//            }


//            for (auto &nfl:node_first_later) /*if (nfl.second.first==1 or nfl.second.second==1)*/ candidates.insert(nfl.first);
//            std::cout<<"--- calling clean_thread_nodepositions ---"<<std::endl;
//            if (clean_thread_nodepositions(thread_node_positions,node_first_later,candidates)) {
//                for (auto nid:nodes) {
//                    auto nv = get_nodeview(nid);
//                    node_first_later[nid] = {0, 0};
//                }
//                for (auto &tnp:thread_node_positions) {
//                    if (thread_nextpos[tnp.first]== tnp.second.size()) thread_nextpos[tnp.first]=-1;
//                    if (thread_nextpos[tnp.first]!=-1) {
//                        ++node_first_later[tnp.second[thread_nextpos[tnp.first]].second].first;
//                        //std::cout<<"incremented first for node "<<tnp.second[thread_nextpos[tnp.first]].second<<" from thread "<<tnp.first<<std::endl;
//                        for (auto i=thread_nextpos[tnp.first]+1;i<tnp.second.size();++i) ++node_first_later[tnp.second[i].second].second;
//                    }
//                }
//                std::cout<<"--- clean_thread_nodepositions finished ---"<<std::endl;
//                continue;
//
//            }
//            std::cout<<"--- clean_thread_nodepositions finished ---"<<std::endl;
            candidates.clear();
            for (auto &nfl:node_first_later) if (nfl.second.first and nfl.second.second==1) candidates.insert(nfl.first);

            sgNodeID_t best_candidate=0;
            for (auto c:candidates) if (best_candidate==0 or node_first_later[best_candidate].first<node_first_later[c].first) best_candidate=c;
            if (best_candidate==0) break;
            bool recovered=false;
            //ATTEMPT 1: does the candidate with most firsts only has one second? then just move up on that thread and make it a candidate
            for (auto &tnp:thread_node_positions) {
                if (thread_nextpos[tnp.first]!=-1 and thread_nextpos[tnp.first]<tnp.second.size()-1 and tnp.second[thread_nextpos[tnp.first]+1].second==best_candidate){
//                    std::cout<<"a single thread seems to have a single node misplaced, correcting"<<std::endl;
                    --node_first_later[tnp.second[thread_nextpos[tnp.first]].second].first;//remove the previous node count as first
                    ++thread_nextpos[tnp.first];
                    ++node_first_later[best_candidate].first;//remove the previous node count as first
                    --node_first_later[best_candidate].second;//remove the previous node count as first
                    recovered=true;
                    break;
                }
            }
            if (recovered) continue;
            //ATTEMPT 2: just find the offending thread and blow it
            uint64_t thread_to_delete=0;
            for (auto &tnp:thread_node_positions) {
                if (thread_nextpos[tnp.first]!=-1 and thread_nextpos[tnp.first]<tnp.second.size()-1 and std::find_if(tnp.second.begin()+thread_nextpos[tnp.first],tnp.second.end(),[&best_candidate](auto &x){return x.second==best_candidate;})!=tnp.second.end()) {
                    --node_first_later[tnp.second[thread_nextpos[tnp.first]].second].first;//remove the first count for the first in the thread
                    for (auto i=thread_nextpos[tnp.first]+1;i<tnp.second.size();++i)
                        --node_first_later[tnp.second[i].second].second;//remove the later count for everyone after, this should include the candidate
                    thread_nextpos[tnp.first]=-1;//make the thread finished

                    break;
                }
            }
            if (recovered) continue;
            //TODO: return the highest % of first vs. later? could be a loopy node that once removed unlocks the order
            break;
        }

        for (auto &tn:thread_nextpos){
            if (tn.second==-1) continue;
            if (thread_node_positions[tn.first][tn.second].second==nid){
                ++tn.second;
                if (tn.second==thread_node_positions[tn.first].size()) tn.second=-1;
                else {
                    ++node_first_later[thread_node_positions[tn.first][tn.second].second].first;
                    --node_first_later[thread_node_positions[tn.first][tn.second].second].second;
                }
            }
        }
        sorted_nodes.emplace_back(nid);
        node_first_later.erase(nid);//node is used!
    }
    return sorted_nodes;
}

int64_t place_node(const std::map<sgNodeID_t, int64_t> &node_positions, const std::map<sgNodeID_t, std::vector<std::pair<sgNodeID_t,int64_t>>> & node_distances, sgNodeID_t nid){
    if (node_distances.count(nid)==0) return INT64_MIN; //unplaced
    std::vector<int64_t> positions;
    for (const auto & nd:node_distances.at(nid)) {
        if (node_positions.count(nd.first))
            positions.emplace_back(node_positions.at(nd.first)+nd.second);
    }
    if (positions.empty()) return INT64_MIN;
    std::sort(positions.begin(),positions.end());
    return positions[positions.size()/2]; //Poor man's median hits back
};

void update_npcomplete(std::map<sgNodeID_t, std::pair<bool,bool>> &np_complete,const std::map<sgNodeID_t, int64_t> &node_positions, const std::map<sgNodeID_t, std::vector<std::pair<sgNodeID_t,int64_t>>> & node_distances, const std::set<sgNodeID_t> &to_place){
    for (auto nid:to_place){
        if (node_distances.count(nid)==0) {
            np_complete[nid]={false,false};
            continue;
        }
        bool prevs=true,nexts=true;
        for (auto nd:node_distances.at(nid)) { //TODO: check prevs/nexts only if still true? saves lookups
            if (node_positions.count(nd.first)==0){
                if (nd.second>0) nexts= false;
                else prevs=false;
            }
        }
        np_complete[nid]={prevs,nexts};
    }
}

sgNodeID_t most_connected_node(const std::map<sgNodeID_t, int64_t> &node_positions, const std::map<sgNodeID_t, std::vector<std::pair<sgNodeID_t,int64_t>>> & node_distances, const std::set<sgNodeID_t> &to_place){
    int most_connected=0;
    sgNodeID_t best_nid=0;
    for (auto nid:to_place){
        if (node_positions.count(nid) or node_distances.count(nid)==0) continue;
        int prevs=0,nexts=0;
        for (auto nd:node_distances.at(nid)) {
            if (node_positions.count(nd.first)){
                if (nd.second>0) ++nexts;
                else ++prevs;
            }
        }
        if (nexts+prevs>most_connected){
            best_nid=nid;
            most_connected=nexts+prevs;
        }
    }
    return best_nid;
}

std::map<sgNodeID_t, std::vector<std::pair<sgNodeID_t,int64_t>>> tnp_to_distances (const std::map<uint64_t, std::vector<std::pair<int64_t, sgNodeID_t>>> &thread_nodepositions,const std::set<sgNodeID_t> &nodeset){
    std::map<sgNodeID_t, std::vector<std::pair<sgNodeID_t,int64_t>>> distances;

    for (const auto &tnp:thread_nodepositions) {
        sgNodeID_t last_node=0;
        int64_t last_pos=0;
        for (auto p:tnp.second){
            if (last_node!=0) {
                if (nodeset.count(p.second)) distances[p.second].emplace_back(last_node,p.first-last_pos);
                if (nodeset.count(last_node)) distances[last_node].emplace_back(p.second,last_pos-p.first);
            }
            last_node=p.second;
            last_pos=p.first;
        }
    }
    return distances;
}

//TODo: this creates tnps for all the places nodes too, which is not needed.
std::vector<std::pair<sgNodeID_t, int64_t>> ReadThreadsGraph::place_nodes(
        const std::vector<std::pair<sgNodeID_t, int64_t>> &placed_nodes, const std::vector<sgNodeID_t> &nodes,
        bool verbose) const {

    if (verbose) std::cout<<"Place nodes starting"<<std::endl;

    //STEP 1: populate node distances for each node
    std::set<sgNodeID_t> nodeset(nodes.begin(),nodes.end());
    auto full_nodeset=nodeset;
    for (auto &pn:placed_nodes) full_nodeset.insert(pn.first);
    if (verbose) std::cout<<"Nodeset has "<<nodeset.size()<<" nodes, full nodeset has "<<full_nodeset.size()<<std::endl;
    auto tnp=make_thread_nodepositions(full_nodeset);
    if (verbose) std::cout<<"Created thread nodepositions from "<<tnp.size()<<" threads"<<std::endl;
    auto node_distances=tnp_to_distances(tnp,nodeset);
    if (verbose) std::cout<<"Created distances for "<<node_distances.size()<<" nodes"<<std::endl;


    //STEP 2: place each node in the median of its placed distances.
    std::vector<std::map<sgNodeID_t, int64_t>> v_node_positions; //All possible starting conditions

    if (not placed_nodes.empty()) { //just one starting position given as parameter
        v_node_positions.emplace_back();
        for (const auto &np:placed_nodes) v_node_positions.back()[np.first] = np.second;
    }
    else { //for every node that has no prevs, use it as starting position in 0
        for (auto &nds:node_distances) {
            bool noprevs=true;
            for (auto nd:nds.second) {
                if (nd.second>0) {
                    noprevs=false;
                    break;
                }
            }
            if (noprevs) {
                v_node_positions.emplace_back();
                v_node_positions.back()[nds.first]=0;
                if (verbose) std::cout<<"creating a starting position with "<<nds.first<<" at 0"<<std::endl;
            }
        }
    }
    int64_t max_placed_nodes=0,best_solution=-1;

    for (int64_t npix=0;npix<v_node_positions.size();++npix){
        auto &node_positions=v_node_positions[npix];
        std::set<sgNodeID_t> to_place(nodes.begin(),nodes.end());
        std::map<sgNodeID_t, std::pair<bool,bool>> np_complete;
        if (verbose) std::cout<<"Entering placing loop"<<std::endl;
        while (not to_place.empty()){
            update_npcomplete(np_complete,node_positions,node_distances,to_place);
            std::set<sgNodeID_t> placed;
            //first place with both prev and nexts full; if none, with prevs full; if none, with nexts full
            for (std::pair<bool,bool> cond:{std::make_pair(true,true),{true,false},{false,true}}) {

                for (auto nid:to_place) {
                    if (np_complete[nid] == cond) {
                        auto p = place_node(node_positions, node_distances, nid);
                        if (p != INT64_MIN) {
                            node_positions[nid] = p;
                            placed.insert(nid);
                        }
                    }
                }
                if (not placed.empty()) break;
            }


            //if none placed, abort (for now), we can change np_complete to count unsatisfied and place the smallest
            if (placed.empty()) {
                //find node with most connections, place that one
                //std::cout<<"Aborting after placing "<<node_positions.size()<<" nodes, with "<<to_place.size()<<" nodes still to place, but no candidates"<<std::endl;
                auto nid=most_connected_node(node_positions,node_distances,to_place);
                auto p = place_node(node_positions, node_distances, nid);
                if (p != INT64_MIN) {
                    node_positions[nid] = p;
                    placed.insert(nid);
                }
                else {
                    if (verbose) std::cout<<"Aborting after placing "<<node_positions.size()<<" nodes, with "<<to_place.size()<<" nodes still to place, but no candidates"<<std::endl;
                    break;
                }
            }
            for (auto nid:placed) to_place.erase(nid);
        }
        if (node_positions.size()>max_placed_nodes) {
            max_placed_nodes=node_positions.size();
            best_solution=npix;
        }
    }
    if (best_solution==-1) return {};
    auto &node_positions=v_node_positions[best_solution];
    std::vector<std::pair<sgNodeID_t,int64_t>> full_placements;
    full_placements.reserve(node_positions.size());
    for (auto &np:node_positions) full_placements.emplace_back(np.first,np.second);
    std::sort(full_placements.begin(),full_placements.end(),[](auto &a,auto &b){return a.second<b.second;});
    if (placed_nodes.empty() and full_placements[0].second!=0){
        auto offset=full_placements[0].second;
        for (auto &np:full_placements) np.second-=offset;
    }
    return full_placements;
}