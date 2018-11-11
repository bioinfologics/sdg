//
// Created by Bernardo Clavijo (EI) on 25/05/2018.
//

#ifndef BSG_LINKAGEDIGRAPH_HPP
#define BSG_LINKAGEDIGRAPH_HPP

#include <sglib/graph/SequenceGraph.hpp>

/**
 * This is a description of the graph creating new linkage, the nodes are referenced from the original graph.
 *
 * The LinkageDIgraph has links describing structure in the graph but not the nodes, the nodes are references from sg original graph,
 * the linkage digraph references to the nodes of a sequence digraph and only contains the links.
 * It's usually used to encode solutions from algorithms.
 *
 * Each link is represented on each node it appears as -A,B on A and -B,A on B
 */
class LinkageDiGraph {
public:
    LinkageDiGraph(SequenceGraph & _sg): sg(_sg){};

    /** @brief Adds a link between source and destination in the links collection.
     * Each link is added from both ends in the collection (see links vector)
     * The link is directed
     *
     * @param source Source node + or - end
     * @param dest Destination node + or - end
     * @param d Distance (overlap if it's negative)
     */
    void add_link( sgNodeID_t source, sgNodeID_t dest, int32_t d);

    /**
     * Given a LDG adds all the connection of that LDG to the current LDG links collection.
     *
     * @param other LDG instance
     */
    void add_links(const LinkageDiGraph &other);

    /** @brief If the link between source and dest exists removes the link from the link collection
     *
     * @param source
     * @param dest
     */
    void remove_link(sgNodeID_t source, sgNodeID_t dest);

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

    /** @brief Given a nodeID returns a set of all fw nodes present up to radius (# jump) away
     * TODO: merge with DFS in sequence graph (?)
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
     * @param nodeID n1
     * @param nodeID n2
     * @return true is the nodes are connected false otherwise
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
     * Dumps the LinkageDiGraph (links) to <filename>
     * @param filename
     */
    void dump_to_text(std::string filename);

    /**
     * Loads the LinkageDiGraph (links) from <filename>
     * @param filename
     */
    void load_from_text(std::string filename);

    SequenceGraph & sg;

    /**
     * Collection of links for each node of the graph
     * links[node] = [link, link....]
     */
    std::vector<std::vector<Link>> links;

};
#endif //BSG_LINKAGEDIGRAPH_HPP
