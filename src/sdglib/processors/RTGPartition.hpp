//
// Created by Bernardo Clavijo (EI) on 06/08/2021.
//

#pragma once
#include <sdglib/graph/ReadThreadsGraph.hpp>
#include <sdglib/processors/HappySorter.hpp>

class RTGPartition {
public:
    RTGPartition(const ReadThreadsGraph &rtg,int min_node_threads=2,float node_min_percentage=.6, int thread_p=10, int thread_q=12, int min_thread_nodes=-1);


    int64_t get_node_class(sgNodeID_t nid);
    int64_t compute_node_class(sgNodeID_t nid);
    bool switch_node_class(sgNodeID_t nid, int64_t c);

    int64_t get_thread_class(int64_t tid);
    int64_t compute_thread_class(int64_t tid, int distance_to_end= 5);
    bool switch_thread_class(int64_t tid, int64_t c);

    std::vector<int64_t> find_unclassified_threads(int min_nodes=50, float max_classified_nodes_perc=.05);

    int64_t new_class_from_thread(int64_t tid);

    LocalOrder order_from_class(int64_t cid);

    bool supported_thread(int64_t tid,int min_support=2);

    void classify_all_threads(int min_nodes=50, float max_classified_nodes_perc=.05);

    uint64_t propagate(uint64_t steps=UINT64_MAX, bool verbose=false);

    std::unordered_map<std::pair<int64_t,int64_t>,std::vector<int64_t>> find_class_bridges(int p, int q);

    void reset(int min_node_threads=-1,float node_min_percentage=-1, int thread_p=-1, int thread_q=-1);

    const ReadThreadsGraph &rtg;
    int min_node_threads;
    float node_min_percentage;
    int thread_p;
    int thread_q;
    int min_thread_nodes;
    std::unordered_map<sgNodeID_t,std::vector<int64_t>> node_threads;
    std::unordered_map<int64_t,std::vector<sgNodeID_t>> thread_nodes;
    std::unordered_map<sgNodeID_t,int64_t> node_class;
    std::unordered_map<int64_t ,int64_t> thread_class;
    std::unordered_set<sgNodeID_t> nodes_to_evaluate;
    std::unordered_set<int64_t> threads_to_evaluate;

private:
    bool thread_available(uint64_t tid, int min_nodes=50, float max_classified_nodes_perc=.05);
};


