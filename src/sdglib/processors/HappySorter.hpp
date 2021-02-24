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

    void reverse();
    float thread_happiness(int64_t tid,int min_nodes=-1) const;
    float node_happiness(sgNodeID_t,bool prev=true,bool next=false, int min_threads=-1) const; //TODO: have an open-thread happiness
    void recruit_all_happy_threads(float min_happiness=-1, int min_nodes=-1);
    void close_internal_threads(int order_end=20,int thread_end=0);
    std::unordered_set<sgNodeID_t> find_fw_candidates(float min_happiness=-1, int min_threads=-1, int end_size=-1) const;
    std::unordered_set<sgNodeID_t> find_bw_candidates(float min_happiness=-1, int min_threads=-1, int end_size=-1) const;
    //This looks for candidates that are both fw for some nodes and bw for others
    std::unordered_set<sgNodeID_t> find_internal_candidates(float min_happiness=-1, int min_threads=-1, int32_t first=1, int32_t last=INT32_MAX) const;
    void start_from_node(sgNodeID_t nid, int min_links=4);
    //this works similarly to the rtg one, but includes all nodes form the order in the threads too
    std::map<int64_t,std::vector<std::pair<int64_t,sgNodeID_t>>> make_thread_nodepositions(const std::unordered_set<sgNodeID_t> &nodes) const;
    std::vector<std::pair<sgNodeID_t, int64_t>> place_nodes( const std::unordered_set<sgNodeID_t> &nodes, bool verbose) const;
    bool  add_placed_nodes( const std::vector<std::pair<sgNodeID_t, int64_t>> &placed_nodes, bool update_current=true);
    bool grow_fw(int min_threads, bool verbose=true);
    bool grow(int min_threads=-1, float min_happiness=-1, bool fw=true, bool bw=true, bool internal=true);
    bool grow_loop(int min_threads=-1, float min_happiness=-1, int64_t steps=INT64_MAX);

    float min_thread_happiness;
    int min_thread_nodes;
    float min_node_happiness;
    int min_node_threads;
    int order_end_size;

    const ReadThreadsGraph & rtg;

    LocalOrder order;

    std::unordered_map<sgNodeID_t,int64_t> node_coordinates;

    std::unordered_set<int64_t> threads;
    std::unordered_set<int64_t> fw_open_threads;
    std::unordered_set<int64_t> bw_open_threads;

};

class HappySorterRunner {

};
