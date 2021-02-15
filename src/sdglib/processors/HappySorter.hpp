//
// Created by Bernardo Clavijo (EI) on 12/02/2021.
//

#pragma once

#include <sdglib/graph/ReadThreadsGraph.hpp>
#include <sdglib/processors/ThreadedGraphSorter.h>

class HappySorter {
public:
    HappySorter(const ReadThreadsGraph & _rtg, float _min_thread_happiness=.7, float _min_node_happiness=.7,
                int _min_thread_nodes=3, int _min_node_threads=2, int _order_end_size=20):
            rtg(_rtg), min_thread_happiness(_min_thread_happiness), min_node_happiness(_min_node_happiness),
            min_thread_nodes(_min_thread_nodes), min_node_threads(_min_node_threads), order_end_size(_order_end_size){};

    //TODO: recruit all happy intermediate nodes
    //TODO: start from node
    //TODO: grow fw
    //TODO: reverse
    //TODO: fill-in
    //TODO: close threads

    float thread_happiness(int64_t tid,int min_nodes=-1) const;
    float node_happiness(sgNodeID_t,bool prev=true,bool next=false, int min_threads=-1) const;
    void recruit_all_happy_threads(float min_happiness=-1, int min_nodes=-1);
    std::unordered_set<sgNodeID_t> find_fw_candidates(float min_happiness=-1, int min_threads=-1, int end_size=-1) const;
    //This looks for candidates that are both fw for some nodes and bw for others
    std::unordered_set<sgNodeID_t> find_internal_candidates(float min_happiness=-1, int min_threads=-1) const;

    float min_thread_happiness;
    int min_thread_nodes;
    float min_node_happiness;
    int min_node_threads;
    int order_end_size;

    const ReadThreadsGraph & rtg;

    LocalOrder order;

    std::unordered_set<int64_t> threads;
    std::unordered_set<int64_t> fw_open_threads;
    std::unordered_set<int64_t> bw_open_threads;

};
