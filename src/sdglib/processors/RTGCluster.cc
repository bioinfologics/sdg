//
// Created by Bernardo Clavijo (EI) on 04/08/2021.
//

#include "RTGCluster.hpp"

RTGCluster::RTGCluster(const ReadThreadsGraph &_rtg, int _p, int _q, int _min_node_threads, float _min_node_happiness):
    rtg(_rtg), p(_p), q(_q), min_node_threads(_min_node_threads), min_node_happiness(_min_node_happiness) {
    int node_hap_perc=_min_node_happiness*100;
    for (const auto &nv:rtg.get_all_nodeviews(false,false)) {
        node_total_threads[nv.node_id()] = node_total_threads[-nv.node_id()] = rtg.node_threads(nv.node_id()).size();
        int ht=node_hap_perc*node_total_threads[nv.node_id()]/100;
        node_happiness_threads[nv.node_id()] = node_happiness_threads[-nv.node_id()] = std::max(ht,min_node_threads);
    }
}

bool RTGCluster::is_node_happy(sgNodeID_t nid) {
    return node_threads[nid]>=node_happiness_threads[nid];
}

std::vector<sgNodeID_t> RTGCluster::new_happy_nodes() {
    std::vector<sgNodeID_t> hn;
    for (const auto &nt:node_threads) {
        if (nodes.count(nt.first) == 0 and is_node_happy(nt.first)) hn.emplace_back(nt.first);
    }
    return hn;
}

void RTGCluster::add_node(sgNodeID_t nid) {
    if (nodes.count(nid) ==0) {
        nodes.insert(nid);
        for (const auto &tid : rtg.node_threads(nid, true))
            ++thread_nodes[tid];
    }
}

bool RTGCluster::is_thread_happy(int64_t tid) {
    if (thread_nodes[tid]<p) return false;
    const auto thread=rtg.get_thread(tid);
    //This uses a "circular queue" over a linear vec<bool>
    std::vector<bool> node_in_cluster(thread.size()+q,false);
    auto bp=node_in_cluster.begin();
    auto fp=node_in_cluster.begin()+q;
    int current_used=0;
    for (auto &np:thread){
        if (*bp) --current_used;
        if (nodes.count(np.node)){
            ++current_used;
            if (current_used>=p) return true;
            *fp=true;
        }
        ++fp,++bp;
    }
    return false;
}

std::vector<int64_t> RTGCluster::new_happy_threads() {
    std::vector<int64_t> ht;
    for (const auto &tn:thread_nodes) {
        if (threads.count(tn.first) == 0 and is_thread_happy(tn.first)) ht.emplace_back(tn.first);
    }
    return ht;
}

void RTGCluster::add_thread(int64_t tid) {
    if (threads.count(tid)==0) {
        threads.insert(tid);
        for (auto const &np:rtg.get_thread(tid)) {
            ++node_threads[np.node];
        }
    }
}

bool RTGCluster::grow(uint64_t steps) {
    auto start_count=nodes.size()+threads.size();
    if (steps) for (auto &nid:new_happy_nodes()) add_node(nid);
    uint64_t last_count=0;
    while (steps-- and last_count<nodes.size()+threads.size()) {
        last_count=nodes.size()+threads.size();
        for (auto &tid:new_happy_threads()) add_thread(tid);
        for (auto &nid:new_happy_nodes()) add_node(nid);
    }
    return nodes.size()+threads.size()>start_count;
}
