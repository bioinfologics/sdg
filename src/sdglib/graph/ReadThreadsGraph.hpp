//
// Created by Bernardo Clavijo (EI) on 30/11/2020.
//

#pragma once

#include <sdglib/mappers/LongReadsRecruiter.hpp>
#include <sdglib/graph/DistanceGraph.hpp>


/**
 * The ReadThreadsGraph is a DistanceGraph where the links form read threads, it has added conditions and functions to
 * work directly with threads, as to make it more efficient
 */
class ReadThreadsGraph : public DistanceGraph {
public:
    explicit ReadThreadsGraph(SequenceDistanceGraph & sdg): DistanceGraph(sdg,"SDG") {};

    //DistanceGraph subgraph_from_node(sgNodeID_t nid, uint32_t max_distance, int min_links);
    bool add_thread(int64_t thread_id,const std::vector<NodePosition> & node_positions, bool remove_duplicated, int min_thread_nodes);
    bool remove_thread(int64_t thread_id);
    // std::vector<std::pair<int64_t,sgNodeID_t>> sort_graph();
    //bool pop_node(sgNodeID_t node_id,int64_t thread_id);
    //bool pop_node_from_all(sgNodeID_t node_id);
    //std::vector<NodePosition> get_thread(int64_t thread_id);
    //bool split_thread_at(int64_t thread_id, int lidx); FUTURE

    std::unordered_map<int64_t,std::pair<sgNodeID_t,sgNodeID_t>> thread_ends;
    //TODO: maybe save the exact positions of nodes in threads to directly compute distances between any two?
};