//
// Created by Bernardo Clavijo (EI) on 30/11/2020.
//

#pragma once

#include <sdglib/graph/DistanceGraph.hpp>
#include <sdglib/mappers/LongReadsRecruiter.hpp>
#include <sdglib/views/NodeView.hpp>


class ThreadInfo{
public:
    sgNodeID_t start=0;
    sgNodeID_t end=0;
    uint16_t link_count=0;
};

/**
 * The ReadThreadsGraph is a DistanceGraph where the links form read threads, it has added conditions and functions to
 * work directly with threads, as to make it more efficient
 */
class ReadThreadsGraph : public DistanceGraph {
public:
    explicit ReadThreadsGraph(SequenceDistanceGraph & sdg,const std::string &_name="RTG"): DistanceGraph(sdg,_name) {};

    void dump(std::string filename);
    void load(std::string filename);

    //ReadThreadsGraph subgraph_from_node(sgNodeID_t nid, uint32_t max_distance, int min_links);
    bool add_thread(int64_t thread_id,const std::vector<NodePosition> & node_positions, bool remove_duplicated=true, int min_thread_nodes=2);
    bool remove_thread(int64_t thread_id);

    NodeView thread_start_nodeview(int64_t thread_id);
    NodeView thread_end_nodeview(int64_t thread_id);
    LinkView next_in_thread(sgNodeID_t nid, int64_t thread_id,int64_t link_index=-1);
    LinkView prev_in_thread(sgNodeID_t nid, int64_t thread_id,int64_t link_index=1);
    std::vector<sgNodeID_t> all_nids_fw_in_thread(sgNodeID_t nid, int64_t thread_id);
    ReadThreadsGraph local_graph(sgNodeID_t nid,int64_t distance,uint16_t min_links);

    std::unordered_set<uint64_t> node_threads(sgNodeID_t nid);
    // std::vector<std::pair<int64_t,sgNodeID_t>> sort_graph();
    bool pop_node(sgNodeID_t node_id,int64_t thread_id);
    bool pop_nodes(std::vector<sgNodeID_t> node_ids, int64_t thread_id);
    bool pop_node_from_all(sgNodeID_t node_id);
    std::vector<NodePosition> get_thread(int64_t thread_id);
    bool flip_thread(int64_t thread_id);
    std::unordered_map<uint64_t,std::set<sgNodeID_t>> thread_nodesets();
    //bool split_thread_at(int64_t thread_id, int lidx); FUTURE
    std::unordered_map<sgNodeID_t,uint64_t> node_thread_neighbours(sgNodeID_t nid);
    int clean_node(sgNodeID_t node_id, int min_supported=4, int min_support=1);
    std::vector<std::pair<uint64_t,sgNodeID_t>> clean_repeat_nodes_popping_list(int max_threads=200);
    std::vector<std::pair<uint64_t,sgNodeID_t>> clean_all_nodes_popping_list(int min_supported=4, int min_support=1);
    std::vector<std::pair<uint64_t,sgNodeID_t>> clean_all_nodes_by_thread_clustering_popping_list(int min_shared=4, float max_second_perc=.1);
    void apply_popping_list(const std::vector<std::pair<uint64_t,sgNodeID_t>> &popping_list);
    std::unordered_map<int64_t,ThreadInfo> thread_info;
    //TODO: maybe save the exact positions of nodes in threads to directly compute distances between any two?
};