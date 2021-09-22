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

    int64_t class_from_sorter(HappySorter & sorter,int64_t cid=0);

    HappySorter sorter_from_class(int64_t cid, float _min_node_happiness=.7, int _min_node_threads=3, int _order_end_size=30);

    void classify_all_threads(int min_nodes=50, float max_classified_nodes_perc=.05);

    uint64_t propagate(uint64_t steps=UINT64_MAX, bool verbose=false);

    uint64_t thread_propagate(uint64_t steps=UINT64_MAX, float vote_perc=.1, int max_noise=3, bool verbose=false);

    std::map<uint64_t, uint64_t> thread_class_votes(uint64_t tid, std::set<ThreadOverlapType> ovltypes={});

    uint64_t thread_propagate2(uint64_t steps=UINT64_MAX, int min_side_votes=3, float side_vote_perc=.75, int min_contain_votes=5, float contain_vote_perc=.90,  bool reclassify= false, bool verbose=false);

    std::unordered_map<std::pair<int64_t,int64_t>,std::vector<int64_t>> find_class_bridges(int p, int q);
    //find_class_bridges -> threads on class 0, where %threads on nodes switches from start to end between two classes
        //could also be done will all thread's p/q, then comparing all p/q's
    //is_class_mixed: -> can use lines to check? also, connectivity/subclasses

    void compute_thread_intersections(int min_threads, int max_threads);

    void compute_thread_neighbours(int min_shared=10);

    void compute_thread_neighbours_p_q(int p=10, int q=10);

    std::vector<int64_t> get_thread_neighbours(int64_t tid) const;

    std::vector<ThreadOverlapType> get_thread_neighbours_types(int64_t tid) const;

    void reset(int min_node_threads=-1,float node_min_percentage=-1, int thread_p=-1, int thread_q=-1);

    int64_t get_thread_intersection(int64_t tid1, int64_t tid2) const;

    std::vector<int> thread_shared_detail(int64_t tid1, int64_t tid2) const;

    void classify_neighbours(int skip_nodes=10, bool discard_invalid=true);

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
    std::unordered_map<int64_t,std::vector<ThreadOverlapType>> thread_neighbours_types;

private:
    bool thread_available(uint64_t tid, int min_nodes=50, float max_classified_nodes_perc=.05);
};


