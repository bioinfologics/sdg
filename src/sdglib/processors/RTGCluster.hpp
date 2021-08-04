//
// Created by Bernardo Clavijo (EI) on 04/08/2021.
//

#pragma once
#include <sdglib/graph/ReadThreadsGraph.hpp>

class RTGCluster {
public:
    RTGCluster(const ReadThreadsGraph & _rtg, int _p, int _q, int _min_node_threads, float _min_node_happiness);
    bool is_node_happy(sgNodeID_t nid);
    std::vector<sgNodeID_t> new_happy_nodes();
    void add_node(sgNodeID_t nid);

    bool is_thread_happy(int64_t tid);
    std::vector<int64_t> new_happy_threads();
    void add_thread(int64_t tid);

    bool grow(uint64_t steps=UINT64_MAX);

    const ReadThreadsGraph & rtg;
    int p,q,min_node_threads;
    float min_node_happiness;

    std::unordered_set<sgNodeID_t> nodes;
    std::unordered_set<int64_t> threads;
    std::unordered_map<sgNodeID_t,uint64_t> node_total_threads;
    std::unordered_map<sgNodeID_t,uint64_t> node_happiness_threads;
    std::unordered_map<sgNodeID_t,uint64_t> node_threads;
    std::unordered_map<int64_t,uint64_t> thread_nodes;
};


