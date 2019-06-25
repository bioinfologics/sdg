//
// Created by Bernardo Clavijo (EI) on 25/05/2018.
//

#ifndef BSG_DISTANCEGRAPH_HPP
#define BSG_DISTANCEGRAPH_HPP
#include <vector>
#include <algorithm>
#include <string>
#include <map>
#include <unordered_map>
#include <set>
#include <iostream>
#include <array>
#include <unordered_set>
#include <iosfwd>
#include <sdglib/types/KmerTypes.hpp>
#include <sdglib/types/GenericTypes.hpp>
#include <sdglib/graph/SequenceSubGraph.hpp>
#include <sdglib/graph/SequenceGraphPath.hpp>
#include <sdglib/logger/OutputLog.hpp>

class SequenceDistanceGraph;//fwd declaration (to break circular dependence)
/**
 * This is a description of the graph creating new linkage, the nodes are referenced from the original graph.
 *
 * The LinkageDIgraph has links describing structure in the graph but not the nodes, the nodes are references from sg original graph,
 * the linkage digraph references to the nodes of a sequence digraph and only contains the links.
 * It's usually used to encode solutions from algorithms.
 *
 * Each link is represented on each node it appears as -A,B on A and -B,A on B
 */
class DistanceGraph {
public:
    explicit DistanceGraph(SequenceDistanceGraph & _sdg): sdg(_sdg){};
    DistanceGraph(SequenceDistanceGraph& sdg, std::ifstream &input_file);
    DistanceGraph(SequenceDistanceGraph & _sdg, const std::string& name) : name(name), sdg(_sdg){}

    /** @brief Adds a link between source and destination in the links collection.
     * Each link is added from both ends in the collection (see links vector)
     * The link is directed
     *
     * @param source Source node + or - end
     * @param dest Destination node + or - end
     * @param d Distance (overlap if it's negative)
     * @param read_id id of the read supporting the link (experimental)
     */
    void add_link( sgNodeID_t source, sgNodeID_t dest, int32_t d, uint64_t read_id=0);

    /**
     * Given a LDG adds all the connection of that LDG to the current LDG links collection.
     *
     * @param other LDG instance
     */
    void copy_links(const DistanceGraph &other);

    /** @brief If the link between source and dest exists removes the link from the link collection
     * TODO: this is broken in a Multi-Distance Graph
     * @param source
     * @param dest
     */
    bool remove_link(sgNodeID_t source, sgNodeID_t dest);

    /** @brief Removes all the links in the collection from and to a given nodeID
     * TODO: check this one
     * @param node
     */
    void disconnect_node(sgNodeID_t node);

    /**
     * Given a nodeID returns a vector of links forward from that node
     * if there is no fw links returns an empty vector
     *
     * @param n nodeID
     * @return vector of links forward of the provided node
     */
    std::vector<Link> get_fw_links( sgNodeID_t n) const;

    /**
     * Given a nodeID returns a vector of links backwards from that node
     * if there is no bw links returns an empty vector
     *
     * @param n nodeID
     * @return vector of links backwards of the provided node
     */
    std::vector<Link> get_bw_links( sgNodeID_t n) const;

    /** DEPRECATED!!!!!!!!!!
     * TODO: delete
     *
     * @param source
     * @param dest
     * @return
     */
    Link get_link(sgNodeID_t source, sgNodeID_t dest);

    /**
     * Find node IDs for forward nodes of n
     * TODO: sort-uniq (for multi distances)
     * @param n Node to search neighbours for
     * @return
     * Returns a vector of forward node IDs
     */
    std::vector<sgNodeID_t> get_next_nodes(sgNodeID_t n) const;

    /**
     * Find node IDs for backward nodes of n
     * TODO: sort-uniq (for multi_distances)
     * @param n Node to search neighbours for
     * @return
     * Returns a vector of backward node IDs
     */
    std::vector<sgNodeID_t> get_prev_nodes(sgNodeID_t n) const;

    /** @brief Given a nodeID returns a set of all fw nodes present up to radius (# jump) away
     * TODO: merge with other similar functions and deprecate
     *
     * @param n nodeID
     * @param radius number of jumps to consider
     * @return set of nodeIDs
     */
    std::set<sgNodeID_t> fw_reached_nodes(sgNodeID_t n, int radius) const;

    /**
     * Vector of FW distances and neighbours from n, sorted in ascending distance (median if multi-link) order.
     * @param n
     * @param min_links
     * @return
     */
    std::vector<std::pair<int,sgNodeID_t>> fw_neighbours_by_distance( sgNodeID_t n, int min_links) const;

    /** @brief Returns all the nodeIDs with connections in the link collection
     *
     * @return set of nodeIDs
     */
    std::unordered_set<sgNodeID_t>  get_connected_nodes() const;

    /**
     * Explore in a growing radius and returns all possible paths within limits between 2 node ends.
     * @param from
     * @param to
     * @param max_size
     * @param max_nodes
     * @param abort_on_loops
     * @return
     */
    std::vector<SequenceGraphPath> find_all_paths_between(sgNodeID_t from,sgNodeID_t to, int64_t max_size, int max_nodes=20, bool abort_on_loops=true) const;

    /** @brief Removes transitive connections form the link collection
     *
     * @param radius
     */
    void remove_transitive_links(int radius);

    /** @brief Prints to screen connectivity statistics for the link collection
     *
     * Prints the number of nodes with linkage in the following categories
     * 1-1, 1-0, 1-N, N-N, 0-0
     */
    void report_connectivity();
    //void solve();

    /** @brief returns true if n1 and n2 are connected in the link collection and false otherwise
     *
     * @param n1 nodeID n1
     * @param n2 nodeID n2
     * @return true if the nodes are connected false otherwise
     */
    bool are_connected(sgNodeID_t n1, sgNodeID_t n2) const;
    uint32_t link_count(sgNodeID_t n1, sgNodeID_t n2) const;

    /**
     * Get all unitigs from the graph, returns a vector of paths where all nodes are connected 1-1
     *
     * @param min_nodes minimum ammount of nodes in the unitig to be considered
     * @param min_total_size min size of the unitig to be considered
     * @return vetor of paths
     */
    std::vector<std::vector<sgNodeID_t>> get_all_lines(uint16_t min_nodes, uint64_t min_total_size=0) const;

    /**
     * Finds perfect bubble structures in the graph ( prev -> [n1 | n2] -> next )and returns a vector of pairs being
     * each element of the pair an alternative phase of the bubble
     *
     * @param min_size
     * @param max_size
     * @return
     */
    std::vector<std::pair<sgNodeID_t,sgNodeID_t>> find_bubbles(uint32_t min_size,uint32_t max_size) const;
    std::vector<sgNodeID_t> find_tips(uint32_t min_size=0,uint32_t max_size=1000000) const;
    std::vector<sgNodeID_t> find_self_loops(uint32_t min_size=0,uint32_t max_size=1000000, bool include_circles=true) const;

    /**
     * Dumps the LinkageDiGraph (links) to filename
     * @param filename
     */
    void dump_to_text(std::string filename);

    /**
     * Loads the LinkageDiGraph (links) from filename
     * @param filename
     */
    void load_from_text(std::string filename);

    void write_to_gfa1(std::string filename, const std::vector<sgNodeID_t> &selected_nodes={}, const std::vector<double> &depths={});
    void write_to_gfa2(std::string filename, const std::vector<sgNodeID_t> &selected_nodes={}, const std::vector<double> &depths={});

    DistanceGraph& operator=(const DistanceGraph &o);

    void read(std::ifstream &input_file);
    void write(std::ofstream &output_file);

    SequenceDistanceGraph & sdg;

    /**
     * Collection of links for each node of the graph
     * links[node] = [link, link....]
     */
    std::vector<std::vector<Link>> links;

    std::string name;

};
#endif //BSG_DISTANCEGRAPH_HPP
