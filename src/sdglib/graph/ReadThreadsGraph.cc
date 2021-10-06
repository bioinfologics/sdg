//
// Created by Bernardo Clavijo (EI) on 30/11/2020.
//

#include "ReadThreadsGraph.hpp"
#include <sdglib/views/NodeView.hpp>
#include <atomic>

ReadThreadsGraph::ReadThreadsGraph(const ReadThreadsGraph &other): DistanceGraph(other.sdg,other.name){
    links=other.links;
    thread_info=other.thread_info;
}

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
    auto tid=llabs(thread_id);
    for (const auto &l:links[llabs(nid)]) {
        if (l.support.id == tid and l.source==-nid and (link_index==-1 or l.support.index==link_index)) {
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
    thread.reserve(ti.link_count+1);
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

std::vector<sgNodeID_t> ReadThreadsGraph::get_thread_nodes(int64_t thread_id, bool oriented) const{
    std::vector<sgNodeID_t> thread;
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
    thread.reserve(ti.link_count+1);
    sgNodeID_t nid=0;
    auto tid=llabs(thread_id);
    do {
        if (nid==0) nid=s;
        else {
            for (const auto &l:links[llabs(nid)]) {
                if (l.support.index==lc and l.support.id == tid and l.source==-nid) {
                    nid=l.dest;
                    break;
                }
            }
            lc+=link_increment;
        }
        thread.emplace_back(oriented ? nid : llabs(nid));
    } while (lc>=0 and lc<ti.link_count);
    return thread;
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

ThreadOverlapType ReadThreadsGraph::classify_thread_overlap(int64_t tid1, int64_t tid2, int skip_nodes) const {
    using ToT=ThreadOverlapType;
    auto t1n= get_thread_nodes(tid1);
    std::set<sgNodeID_t> t1s(t1n.begin(),t1n.end());
    auto t2n= get_thread_nodes(tid2);
    std::set<sgNodeID_t> t2s(t2n.begin(),t2n.end());
    std::vector<sgNodeID_t> t2rn;
    t2rn.reserve(t2n.size());
    for (auto rit=t2n.crbegin();rit!=t2n.crend();++rit) t2rn.emplace_back(-*rit);
    std::set<sgNodeID_t> t2rs(t2rn.begin(),t2rn.end());
    auto shared=nodeset_intersection_size(t1s,t2s);
    auto rshared=nodeset_intersection_size(t1s,t2rs);
    bool t2reversed=false;
    if (shared==rshared) return ThreadOverlapType::invalid;
    if (shared<rshared) {
        shared=rshared;
        t2n=t2rn;
        t2s=t2rs;
        t2reversed=true;
    }
    auto to_check=shared;
    int last_shared_i1=-1,last_shared_i2=-1;
    bool t1start=false,t2start=false,t1end=false,t2end=false;
    auto i2=0;
    for (auto i1=0;to_check and i1<t1n.size();++i1){
        if (t2s.count(t1n[i1])){
            if (last_shared_i1!=-1 and i1-last_shared_i1>skip_nodes) return ToT::invalid;
            while (i2<t2n.size() and t2n[i2]!=t1n[i1]) ++i2;
            if (i2==t2n.size()) return ToT::invalid; //there was a permutation in node order
            if (last_shared_i2!=-1 and i2-last_shared_i2>skip_nodes) return ToT::invalid;
            --to_check;
            if (i1<skip_nodes) t1start=true;
            if (i1+skip_nodes>t1n.size()) t1end=true;
            if (i2<skip_nodes) t2start=true;
            if (i2+skip_nodes>t2n.size()) t2end=true;
            last_shared_i1=i1;
            last_shared_i2=i2;
        }
    }
    //if we reached here, overlap is valid on both threads, and ends are set correctly
    if ((not t2start and not t1start) or (not t2end and not t1end)) {
        return ToT::invalid;//repeat (i.e. collase in the middle or Y)
    }
    if (t1start and t1end) {
        if (t2start and t2end) return (t2reversed ? ToT::rcomplete : ToT::complete);
        else return (t2reversed ? ToT::t1_in_rt2 : ToT::t1_in_t2);
    }
    if (t2start and t2end) return (t2reversed ? ToT::rt2_in_t1 : ToT::t2_in_t1);
    if (t1start and t2end) return (t2reversed ? ToT::rt2_t1 : ToT::t2_t1);
    if (t2start and t1end) return (t2reversed ? ToT::t1_rt2 : ToT::t1_t2);
    return ToT::invalid;
}

ReadThreadsGraph ReadThreadsGraph::reduced_graph(int min_thread_nodes, int min_node_threads) {
    //TODO: this can probably be speed up with a node-to-threads map and a threads-to-nodes map, and counters, iterating over the map's vectors, then using the map to filter the graph only once
    ReadThreadsGraph rrtg(*this);
    bool nodes_deleted=true;
    bool threads_deleted=true;
    while (nodes_deleted or threads_deleted) {
        nodes_deleted=false;
        threads_deleted=false;
        ReadThreadsGraph new_rrtg(sdg);
        std::unordered_set<sgNodeID_t> nodes_to_keep;
        for (sgNodeID_t nid=0;nid<rrtg.links.size(); ++nid) if (rrtg.links[nid].size()>0 and rrtg.node_threads(nid).size()>=min_node_threads) {
            nodes_to_keep.insert(nid);
            new_rrtg.links[nid].reserve(rrtg.links[nid].size());
        }
        for (auto &ti:rrtg.thread_info) {
            if (ti.second.link_count<=min_thread_nodes) {
                threads_deleted=true;
                continue;
            }
            auto t=rrtg.get_thread(ti.first);
            std::vector<NodePosition> filtered_thread;
            filtered_thread.reserve(t.size());
            for (auto &tp:t) {
                if (nodes_to_keep.count(llabs(tp.node))) filtered_thread.emplace_back(tp);
                else nodes_deleted=true;
            }
            if (filtered_thread.size()>=min_thread_nodes) new_rrtg.add_thread(ti.first,filtered_thread);
            else threads_deleted=true;
        }
        rrtg.links=std::move(new_rrtg.links);
        rrtg.thread_info=std::move(new_rrtg.thread_info);
    }
    return std::move(rrtg);
}

ReadThreadsGraph ReadThreadsGraph::reduced_graph(std::unordered_set<sgNodeID_t> nodes) {
    ReadThreadsGraph rrtg(sdg);
    for (auto &ti:thread_info) {
        auto t=get_thread(ti.first);
        std::vector<NodePosition> filtered_thread;
        filtered_thread.reserve(t.size());
        for (auto &tp:t) {
            if (nodes.count(tp.node) or nodes.count(-tp.node)) filtered_thread.emplace_back(tp);
        }
        if (filtered_thread.size()>=1) rrtg.add_thread(ti.first,filtered_thread);
    }
    return std::move(rrtg);
}

void ReadThreadsGraph::compute_node_proximity(int radius) {
    //TODO: OpenMP is your friend
    std::unordered_set<std::pair<sgNodeID_t,sgNodeID_t>> node_pairs;
    ++radius;
    for (auto ti:thread_info) {
        auto t= get_thread_nodes(ti.first);
        for (auto i=0;i<t.size();++i) {
            for (auto j=i+1;j<t.size() and j<i+radius;++j) {
                auto n1=(llabs(t[i])<llabs(t[j]) ? t[i] : t[j]);
                auto n2=(llabs(t[i])<llabs(t[j]) ? t[j] : t[i]);
                if (n1<0) {
                    n1=-n1;
                    n2=-n2;
                }
                node_pairs.insert({n1,n2});
            }
        }
    }
    node_proximal_nodes.clear();
    for (auto np:node_pairs) {
        auto n1=np.first;
        auto n2=np.second;
        auto s=shared_threads(n1,n2,true);
        node_proximal_nodes[n1].emplace_back(n1 > 0 ? n2:-n2,s);
        node_proximal_nodes[llabs(n2)].emplace_back(n2 > 0 ? n1:-n1,s);

    }
}

std::unordered_map<sgNodeID_t, uint64_t> ReadThreadsGraph::get_proximal_nodes(sgNodeID_t nid, bool oriented) {
    std::unordered_map<sgNodeID_t, uint64_t> s;
    for (auto &nc:node_proximal_nodes[llabs(nid)]) {
        if (oriented) {
            if (nid>0) s[nc.first]=nc.second;
            else s[-nc.first]=nc.second;
        }
        else s[llabs(nc.first)]+=nc.second;
    }
    return s;
}

int64_t ReadThreadsGraph::shared_threads(sgNodeID_t n1, sgNodeID_t n2, bool oriented) {
    return sdglib::intersection_size(node_threads(n1,oriented),node_threads(n2,oriented));
}