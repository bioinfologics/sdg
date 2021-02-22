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
 * work directly with threads, as to make it more efficient.
 *
 * The thread graph starts empty, threads are added node by node according to need from a donor read thread graph,
 * this graph is usually created using the function
 *
 */
class ReadThreadsGraph : public DistanceGraph {
public:
    explicit ReadThreadsGraph(SequenceDistanceGraph & sdg,const std::string &_name="RTG"): DistanceGraph(sdg,_name) {};

    void dump(std::string filename);
    void load(std::string filename);

    //ReadThreadsGraph subgraph_from_node(sgNodeID_t nid, uint32_t max_distance, int min_links);
    /** @brief Add thread to the RTG
     *
     * @param thread_id id of the thread (rid)
     * @param node_positions Vector of NodePosition of the thread
     * @param remove_duplicated if True the thread is deduplicated before is added to the graph
     * @param min_thread_nodes Min number of nodes for the thread to be added
     * @return
     */
    bool add_thread(int64_t thread_id,const std::vector<NodePosition> & node_positions, bool remove_duplicated=true, int min_thread_nodes=2);

    /** @brief Remove a thread from the graph
     *
     * @param thread_id id of the thread to be removed
     * @return
     */
    bool remove_thread(int64_t thread_id);

    /** @brief Get the NodeView at the start of a thread
     *
     * @param thread_id id of the thread
     * @return nodeview of the start of the selected thread
     */
    NodeView thread_start_nodeview(int64_t thread_id) const;

    /** @brief Get the NodeView at the end of the thread
     *
     * @param thread_id id of the thread
     * @return NodeView of the and of the thread
     */
    NodeView thread_end_nodeview(int64_t thread_id) const;

    /** @brief Get the next link from a node in a thread
     *
     * If there is more than one or no links throws an exception
     *
     * @param nid node id to get the next index from
     * @param thread_id thread id to get the links from
     * @param link_index order of the link to get (default -1 gets any order)
     * @return LinkView describing the next link
     */
    LinkView next_in_thread(sgNodeID_t nid, int64_t thread_id,int64_t link_index=-1) const;

    /** @brief Same as next_in_thread but with prev
     *
     * @param nid node id to get the prev index from
     * @param thread_id thread id to get the links from
     * @param link_index order of the link to get (default -1 gets any order)
     * @return LinkView describing the prev link
     */
    LinkView prev_in_thread(sgNodeID_t nid, int64_t thread_id,int64_t link_index=1) const;

    /** @brief From a thread id get all node ids starting from a node.
     *
     * If nid not in thread_id or nid not exactly once in the selected thread, returns an empty vector
     * @param nid node id to start from
     * @param thread_id thread id
     * @return vector with the ids of all nodes fw from nid in thread_id
     */
    std::vector<sgNodeID_t> all_nids_fw_in_thread(sgNodeID_t nid, int64_t thread_id);

    /** @brief Construct the neighbouring thread graph starting from a node
     *
     * @param nid central node
     * @param distance Maximum distance to include a node counting from the central node
     * @param min_links number of links in the thread graph to include a node
     * @return ReadThreadsGraph graph
     */
    ReadThreadsGraph local_graph(sgNodeID_t nid,int64_t distance,uint16_t min_links);

    /** @brief Get all thread ids that go across the node
     *
     * @param nid node id to get the threads from
     * @return set of thread ids
     */
    std::unordered_set<int64_t> node_threads(sgNodeID_t nid,bool oriented=false) const;
    std::unordered_set<std::pair<int64_t,int32_t>> node_threadpositions(sgNodeID_t nid) const;
    // std::vector<std::pair<int64_t,sgNodeID_t>> sort_graph();

    /** @brief Removes the node nid from the thread thread_id in the graph (inplace)
     *
     * @param node_id node to be removed from the thread, the function is not sign sensitive, removes the node in any
     * direction
     * @param thread_id id of the thread to remove the node from
     * @return true is node was removed false otherwise
     */
    bool pop_node(sgNodeID_t node_id,int64_t thread_id);

    /** @brief Removes a group of nodes from a thread in the graph (inplace)
     *
     * @param node_ids vector of node ids to be removed
     * @param thread_id thread_id id of the thread to remove the node from
     * @return true if at least one of the nodes was removed, otherwise false
     */
    bool pop_nodes(std::vector<sgNodeID_t> node_ids, int64_t thread_id);

    /** @brief Pop the node from all the threads in the graph (inplace)
     *
     * The pop_node function was executed for all
     * threads where the node is present).
     * @param node_id Node id to be popped
     * @return false if the node is not present in any thread true otherwise
     */
    bool pop_node_from_all(sgNodeID_t node_id);

    /** @brief Get a thread as an ordered vector of node positions describing the thread
     *
     * @param thread_id thread id
     * @return vector of NodePositions for that specific thread
     */
    std::vector<NodePosition> get_thread(int64_t thread_id) const;

    bool flip_thread(int64_t thread_id);

    /** @brief Creates a map with threads as keys and sets of nodes as values
     *
     * @return thread nodeset map
     */
    std::unordered_map<uint64_t,std::set<sgNodeID_t>> thread_nodesets();
    //bool split_thread_at(int64_t thread_id, int lidx); FUTURE

    /** @brief Get all the reaching nodes from a graph with the threads
     *
     * Get all threads where the node is present, get all nodes in those threads and fill the map
     *
     * @param nid node of interest id
     * @return map where the keys are the neigbouring nodes and the values are the link count for each node in the
     * thread graph
     */
    std::unordered_map<sgNodeID_t,uint64_t> node_thread_neighbours(sgNodeID_t nid, bool oriented=false);

    /** @brief Cleans the node of unsupported connections
     *
     * Removes the node from the threads where the links of the node in the thread are not properly supported.
     * The presence of the node in each thread is tested for support with the rest of the connections (using the
     * node_thread_neighbours links map). If a link in the thread is supported by other threads (>min_support) the link
     * is considered supported.
     *
     * If the placement of a node in a threads has too many unsupported links (>min_supported) the node is considered
     * out of place in a thread and popped.
     *
     * @param node_id node to be cleaned
     * @param min_supported min number of supported links for a node to be kept in a thread
     * @param min_support min number of threads that need to support a link to be considered valid
     * @return number of popped elements (remove of the links is done inplace)
     */
    int clean_node(sgNodeID_t node_id, int min_supported=4, int min_support=1);

    /** @brief Creates a list of repetitive nodes to be removed from the threads graph
     *  A node is here considered repetitive if it's connected more than max_threads times in the threads graph.
     *
     * @param max_threads number of threads to consider a node repetitive
     * @return vector of pairs (thread, nodes) to eliminate from the graph
     */
    std::vector<std::pair<uint64_t,sgNodeID_t>> clean_repeat_nodes_popping_list(int max_threads=200);

    /** @brief Creates a list of unsupported nodes to be removed from the threads graph
     * For the definition of unsupported node look uses the same definition as the clean_node function in this class
     *
     * the cleanup is done for all nodes in the graph.
     *
     * @param min_supported min number of supported links for a node to be kept in a thread
     * @param min_support min number of threads that need to support a link to be considered valid
     * @return vector of pairs (thread, nodes) to eliminate from the graph
     */
    std::vector<std::pair<uint64_t,sgNodeID_t>> clean_all_nodes_popping_list(int min_supported=4, int min_support=1);
    std::vector<std::pair<uint64_t,sgNodeID_t>> clean_all_nodes_by_thread_clustering_popping_list(int min_shared=4, float max_second_perc=.1);

    /** @brief Applies popping lists created using  clean_repeat_nodes_popping_list or clean_all_nodes_popping_list
     *
     * Takes a vector if pairs<tid, nid> and pops all those combinations from the graph
     * @param popping_list list of nodes and threads combinations to pop
     */
    void apply_popping_list(const std::vector<std::pair<uint64_t,sgNodeID_t>> &popping_list);

    bool thread_fw_in_node(int64_t tid,sgNodeID_t nid) const;

    std::map<uint64_t,std::vector<std::pair<int64_t,sgNodeID_t>>> make_thread_nodepositions(const std::set<sgNodeID_t> &nodes) const;

    std::map<sgNodeID_t,std::pair<uint64_t,uint64_t>> make_node_first_later(const std::map<uint64_t,std::vector<std::pair<int64_t,sgNodeID_t>>> &thread_node_positions, const std::map<uint64_t,int64_t> &thread_nextpos={});
    bool clean_thread_nodepositions(std::map<uint64_t,std::vector<std::pair<int64_t,sgNodeID_t>>> &thread_node_positions,
                                                      std::set<sgNodeID_t> nodes_to_review);
    std::vector<sgNodeID_t> order_nodes(const std::vector<sgNodeID_t> nodes, bool write_detailed_log=false) const;

    //places nodes in a linear space already containing placed_nodes, returns the positions for all nodes, including those already placed.
    std::vector<std::pair<sgNodeID_t,int64_t>> place_nodes(const std::vector<std::pair<sgNodeID_t,int64_t>> &placed_nodes, const std::vector<sgNodeID_t> &nodes, bool verbose=false) const;

    /**
     * This map stores the information for all the threads of the graph
     * key: thread id
     * value: ThreadInfo object with the thread information (start, end, link_count)
     */
    std::unordered_map<int64_t,ThreadInfo> thread_info;
    //TODO: maybe save the exact positions of nodes in threads to directly compute distances between any two?
};