//
// Created by Bernardo Clavijo (EI) on 06/08/2021.
//

#pragma once
#include <sdglib/graph/ReadThreadsGraph.hpp>

class RTGClassifier {
public:
    RTGClassifier(const ReadThreadsGraph &rtg,int min_node_threads=2,float node_min_percentage=.6, int thread_p=10, int thread_q=12, int min_thread_nodes=-1);

    int64_t get_node_class(sgNodeID_t nid);
    int64_t compute_node_class(sgNodeID_t nid);
    bool switch_node_class(sgNodeID_t nid, int64_t c);

    int64_t get_thread_class(int64_t tid);
    int64_t compute_thread_class(int64_t tid);
    bool switch_thread_class(int64_t tid, int64_t c);

    uint64_t propagate(uint64_t steps=UINT64_MAX, bool verbose=false);

    std::unordered_map<std::pair<int64_t,int64_t>,std::vector<int64_t>> find_class_bridges(int p, int q);
    //find_class_bridges -> threads on class 0, where %threads on nodes switches from start to end between two classes
        //could also be done will all thread's p/q, then comparing all p/q's
    //is_class_mixed: -> can use lines to check? also, connectivity/subclasses

    void compute_thread_intersections(int min_threads, int max_threads);

    void compute_thread_neighbours(int min_shared=10);

    std::vector<int64_t> get_thread_neighbours(int64_t tid) const;

    void reset(int min_node_threads=-1,float node_min_percentage=-1, int thread_p=-1, int thread_q=-1);

    int64_t get_thread_intersection(int64_t tid1, int64_t tid2) const;

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
    std::unordered_map<std::pair<int64_t, int64_t> ,int64_t> thread_intersections;
    std::unordered_map<int64_t,std::vector<int64_t>> thread_neighbours;
};


