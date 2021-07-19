//
// Created by Bernardo Clavijo (EI) on 12/02/2021.
//

#pragma once

#include <sdglib/graph/ReadThreadsGraph.hpp>

class LocalOrder{
public:
    LocalOrder(){};
    explicit LocalOrder(const std::vector<sgNodeID_t> & nodes);
    std::vector<sgNodeID_t> as_signed_nodes() const;
    std::vector<NodePosition> as_thread(const DistanceGraph &dg) const;
    int64_t get_node_position(sgNodeID_t nid) const;
    int64_t place_node(const ReadThreadsGraph &rtg, sgNodeID_t nid, const std::unordered_map<sgNodeID_t,int64_t> & node_coordinates={}, int max_hops=3, int min_links=3);
    bool  add_placed_nodes( const std::vector<std::pair<sgNodeID_t, int64_t>> &placed_nodes, bool update_current=true);
    LocalOrder reverse() const;
    int64_t thread_order_crosses(const std::vector<NodePosition> & thread) const;
    auto size() const {return node_positions.size();};
    LocalOrder merge(const LocalOrder & other,int max_overhang=4,float min_shared_perc=.5,int min_shared=20,float max_disordered_perc=.02) const;

    void read(std::ifstream &ifs);
    void write(std::ofstream &ofs);

    std::unordered_map<sgNodeID_t , int64_t > node_positions;
    std::unordered_map<sgNodeID_t,int64_t> node_coordinates;
};

class HappySorter {
public:
    //***** Updated functions, part of Jun2021 algorithms
    HappySorter(const ReadThreadsGraph & _rtg, float _min_node_happiness=.7, int _min_node_threads=3, int _order_end_size=30):
            rtg(_rtg), min_node_happiness(_min_node_happiness), min_node_threads(_min_node_threads), order_end_size(_order_end_size){};

    float node_happiness(sgNodeID_t,bool prev=true,bool next=false, int min_threads=-1) const;




    bool grow(int _thread_hits=4, int _end_size=-1, int _node_hits=3, float _min_happiness=.7);

    void start_from_nodelist(std::vector<sgNodeID_t> nodes, int min_links=1);

    void read(std::ifstream &ifs);
    void write(std::ofstream &ofs);

    //***** Old functions, still to review and update
    /**@brief
     * Happines is calculated as the propotion of shared threads between the node and the order
     * */

    void close_internal_threads(int order_end=30,int thread_end=0);
    std::unordered_set<sgNodeID_t> find_fw_candidates(float min_happiness=-1, int min_threads=-1, int end_size=-1) const;
    std::unordered_set<sgNodeID_t> find_bw_candidates(float min_happiness=-1, int min_threads=-1, int end_size=-1) const;
    //This looks for candidates that are both fw for some nodes and bw for others
    std::unordered_set<sgNodeID_t> find_internal_candidates(float min_happiness=-1, int min_threads=-1, int32_t first=1, int32_t last=INT32_MAX) const;
    //this works similarly to the rtg one, but includes all nodes form the order in the threads too
    std::map<int64_t,std::vector<std::pair<int64_t,sgNodeID_t>>> make_thread_nodepositions(const std::unordered_set<sgNodeID_t> &nodes,std::set<int64_t> tids={}) const;
    std::vector<std::pair<sgNodeID_t, int64_t>> place_nodes( const std::unordered_set<sgNodeID_t> &nodes, bool verbose) const;
    std::vector<std::pair<sgNodeID_t, int64_t>> place_nodes_ltr( const std::unordered_set<sgNodeID_t> &nodes, sgNodeID_t first_node, bool verbose) const;

    bool q_grow_loop(int min_threads=-1, float min_happiness=-1, int p=4, int q=5, int64_t steps=INT64_MAX, bool verbose=false);

    void recruit_all_happy_threads_q(int min_nodes, int max_span);
    bool thread_happiness_q(int64_t tid,int min_nodes, int max_span) const;

    int64_t hs_place_node(const std::unordered_map<sgNodeID_t, int64_t> &node_positions, const std::map<sgNodeID_t, std::vector<std::pair<sgNodeID_t,int64_t>>> & node_distances, sgNodeID_t nid) const;
    void hs_update_npcomplete(std::map<sgNodeID_t, std::pair<bool,bool>> &np_complete,const std::unordered_map<sgNodeID_t, int64_t> &node_positions, const std::map<sgNodeID_t, std::vector<std::pair<sgNodeID_t,int64_t>>> & node_distances, const std::unordered_set<sgNodeID_t> &to_place) const;
    sgNodeID_t hs_most_connected_node(const std::unordered_map<sgNodeID_t, int64_t> &node_positions, const std::map<sgNodeID_t, std::vector<std::pair<sgNodeID_t,int64_t>>> & node_distances, const std::unordered_set<sgNodeID_t> &to_place) const;
    std::map<sgNodeID_t, std::vector<std::pair<sgNodeID_t,int64_t>>> hs_tnp_to_distances (const std::map<int64_t, std::vector<std::pair<int64_t, sgNodeID_t>>> &thread_nodepositions,const std::unordered_set<sgNodeID_t> &nodeset) const;

    void run_from_nodelist(std::vector<sgNodeID_t> nodes, int min_threads=-1, float min_happiness=-1, int p=4, int q=5, int64_t steps=INT64_MAX, bool verbose=false);

    bool update_positions(int64_t first=0, int64_t last=-1);

    bool is_mixed(int win=50,float fail=.2);

    float min_node_happiness;
    int min_node_threads;
    int order_end_size;

    const ReadThreadsGraph & rtg;

    LocalOrder order;



    std::unordered_set<int64_t> threads;
    std::unordered_set<int64_t> fw_open_threads;
    std::unordered_set<int64_t> bw_open_threads;

};

class HappySorterRunner {
public:
    HappySorterRunner(const ReadThreadsGraph & _rtg, float _min_thread_happiness=.7, float _min_node_happiness=.7,
    int _min_thread_nodes=3, int _min_node_threads=2, int _order_end_size=20):
    rtg(_rtg), min_thread_happiness(_min_thread_happiness), min_node_happiness(_min_node_happiness),
    min_thread_nodes(_min_thread_nodes), min_node_threads(_min_node_threads), order_end_size(_order_end_size){
        node_sorted.resize(rtg.sdg.nodes.size());
    };
    void run_from_dg_lines(DistanceGraph dg, int min_line_nodes=10, int p=4, int q=6);

    void load(std::string filename);
    void dump(std::string filename);

    float min_thread_happiness;
    int min_thread_nodes;
    float min_node_happiness;
    int min_node_threads;
    int order_end_size;

    const ReadThreadsGraph & rtg;
    std::vector<bool> node_sorted;
    std::unordered_map<sgNodeID_t,LocalOrder> orders;
    std::unordered_map<sgNodeID_t,std::vector<sgNodeID_t>> node_orders;
};
