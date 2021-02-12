//
// Created by Bernardo Clavijo (EI) on 12/02/2021.
//

#include "HappySorter.hpp"

float HappySorter::thread_happiness(int64_t tid,int min_nodes) const {
    if (min_nodes==-1) min_nodes=min_thread_nodes;
    int64_t first_ti=-1,last_ti=-1,first_oi=-1,last_oi=-1,shared_nodes=0;
    auto tnp=rtg.get_thread(tid);
    for (auto i=0;i<tnp.size();++i) {
        auto p=order.get_node_position(tnp[i].node);
        if (p>0){
            if (first_ti==-1){
                first_ti=i;
                first_oi=p;
            }
            last_ti=i;
            last_oi=p;
            ++shared_nodes;
        }
    }
    if (last_ti==-1 or shared_nodes<min_nodes) {
        std::cout<<"Thread happiness for thread "<<tid<<" = 0, only "<<shared_nodes<<"/"<<min_nodes<<" requered shared nodes"<<std::endl;
        return 0;
    }
    std::cout<<"Thread happiness for thread"<<tid<<" = "<<((float) shared_nodes)/(last_ti+1-first_ti)<<std::endl<<" based on "<<shared_nodes<<" shared nodes"<<std::endl;
    return ((float) shared_nodes)/(last_ti+1-first_ti);
}

float HappySorter::node_happiness(sgNodeID_t nid, bool prev, bool next,int min_threads) const {
    if (min_threads==-1) min_threads=min_node_threads;
    uint64_t happy=0,total=0;
    if (nid==0 or rtg.sdg.links.size()-1<llabs(nid)) return 0;
    if (not prev and not next) return 0;
    else if (prev and not next) {
        for (auto const &l:rtg.links[llabs(nid)]) {
            if (l.source==-nid) {
                ++total;
                if (threads.count(rtg.thread_fw_in_node(l.support.id,nid) ? l.support.id : -l.support.id)) ++happy;
            }
        }
    }
    else if (not prev and next) {
        for (auto const &l:rtg.links[llabs(nid)]) {
            if (l.source==nid) {
                ++total;
                if (threads.count(rtg.thread_fw_in_node(l.support.id,nid) ? l.support.id : -l.support.id)) ++happy;
            }
        }
    }
    else {
        for (auto const &tid:rtg.node_threads(nid,true)){
            ++total;
            if (threads.count(tid)) ++happy;
        }

    }
    if (total<min_threads) return 0;
    return ((float) happy)/total;
}

void HappySorter::recruit_all_happy_threads(float min_happiness, int min_nodes) {
    if (min_nodes==-1) min_nodes=min_thread_nodes;
    if (min_happiness==-1) min_happiness=min_thread_happiness;
    for (const auto &np: order.node_positions){
        const auto &nid=np.second>0 ? np.first:-np.first;
        for (auto &tid:rtg.node_threads(nid,true)){
            if (threads.count(tid)==0 and threads.count(-tid)==0 and thread_happiness(tid,min_nodes)>=min_happiness) {
                threads.insert(tid);
                fw_open_threads.insert(tid);
                bw_open_threads.insert(tid);
            }
        }
    }
}


std::unordered_set<sgNodeID_t> HappySorter::find_fw_candidates(float min_happiness, int min_threads, int end_size) const {
    if (min_threads==-1) min_threads=min_node_threads;
    if (min_happiness==-1) min_happiness=min_node_happiness;
    if (end_size==-1) end_size=order_end_size;
    std::unordered_map<sgNodeID_t,int> node_links;
    std::unordered_set<sgNodeID_t> candidates;


    for (auto tid:fw_open_threads){
        auto tnps=rtg.get_thread(tid);
        int64_t last_in_node=-1;
        for (auto i=0;i<tnps.size();++i) {
            if (order.get_node_position(tnps[i].node)>order.size()-end_size) {
                last_in_node=i;
                break;
            }
        }
        for (auto i=last_in_node;i<tnps.size();++i) ++node_links[tnps[i].node];
    }

    for (auto &nl:node_links) if (nl.second>=min_threads /*and order.get_node_position(nl.first)==0*/ and node_happiness(nl.first,true,false,min_threads)>=min_happiness) candidates.insert(nl.first);

    //remove candidates that appear in both directions (should be rare?)
    std::vector<sgNodeID_t> to_delete;
    for (auto c:candidates) if (c>0 and candidates.count(-c)) to_delete.push_back(c);
    for (auto cd:to_delete) {
        candidates.erase(cd);
        candidates.erase(-cd);
    }
    return candidates;
}