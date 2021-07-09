//
// Created by Bernardo Clavijo (EI) on 30/06/2021.
//

#pragma once

#include <sdglib/graph/ReadThreadsGraph.hpp>
#include "HappySorter.hpp"

class TotalSorter {
public:
    //start by choosing a list of anchors and threads to place (i.e. the largest set of anchors with x threads and threads with x anchors
    TotalSorter(ReadThreadsGraph &rtg, int min_thread_length=6, int min_node_threads=3);
    TotalSorter(ReadThreadsGraph &rtg, std::string);

    void prune_rtg();
    //use happysorters to create a set of orders, every node and thread belongs to a number of them, classify equivalent ones
    void run_sorters_from_lines(std::vector<std::vector<sgNodeID_t>> _lines,int min_line_size, int min_order_size=200,float line_occupancy=.7);

    void update_usage();
    void remove_mixed(int win=50,float fail=.2);
    //every node and thread can only belong to an equivalent class
    void compute_sorter_classes();
    //process equivalent classes into single orders
    //done!

    void load(std::string filename);
    void dump(std::string filename);

    int next_sorter=1;
    ReadThreadsGraph &rtg;
    std::unordered_set<sgNodeID_t> nodes;
    std::unordered_set<int64_t> threads;
    std::unordered_map<int64_t,HappySorter> sorters;
    std::unordered_map<int64_t,int64_t> sorter_classes;
    std::unordered_map<sgNodeID_t,std::vector<int64_t>> node_sorters;
    std::unordered_map<int64_t ,std::vector<int64_t>> thread_sorters;

};
