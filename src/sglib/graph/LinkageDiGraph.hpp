//
// Created by Bernardo Clavijo (EI) on 25/05/2018.
//

#ifndef BSG_LINKAGEDIGRAPH_HPP
#define BSG_LINKAGEDIGRAPH_HPP

#include <sglib/graph/SequenceGraph.hpp>

/**
 * This is a description of the graph creating new linkage, is like augmentation of the graph using the nodes but no the nodes,
 * the nodes are referenced from the original graph
 *
 * The difference is that sg has nodes and the sequence and also the links corresponding to that original graph, the linkage digraph
 * references to the nodes of a sequence digraph and only contains the links. it's usually used to encode solutions from algorithms
 * it can be seen as a representation of a solution projected over the original nodes.
 *
 *
 * Contains a graph and links generated between nodes on that graph
 * Each link is represented on each node it appears as -A,B on A and -B,A on B
 */
class LinkageDiGraph {
public:
    LinkageDiGraph(SequenceGraph & _sg): sg(_sg){};

    /** @brief Adds a link between source and destination in the links collection.
     * Each link is added from both ends.
     * The link is directed
     *
     * @param source Source node + or - end
     * @param dest Destination node + or - end
     * @param d Distance (overlap if it's negative)
     */
    void add_link( sgNodeID_t source, sgNodeID_t dest, int32_t d);

    /**
     * Given a LDG adds the connection of that LDG to the current LDG links collection.
     *
     * @param other LDG instance
     */
    void add_links(const LinkageDiGraph &other);

    /**
     * Removes a link from the link conection
     *
     * @param source
     * @param dest
     */
    void remove_link(sgNodeID_t source, sgNodeID_t dest);

    /** @brief Removes all the links in the collection from and to a given node
     * TODO: check this one
     * @param node
     */
    void disconnect_node(sgNodeID_t node);

    /**
     * Given a node id returns a vector of links forward from that node
     * if there is no fw links returns an empty vector
     *
     * @param n node id
     * @return vector of links forward of the provided node
     */
    std::vector<Link> get_fw_links( sgNodeID_t n) const;

    /**
     * Given a node id returns a vector of links backwards from that node
     * if there is no bw links returns an empty vector
     *
     * @param n
     * @return
     */
    std::vector<Link> get_bw_links( sgNodeID_t n) const;

    /** @brief Given a node id returns a set of all fw nodes present up to radius (# jump) away
     * TODO: unfy with DFS in sequence graph (?)
     * @param n node id
     * @param radius number of jumps to consider
     * @return set of node ids
     */
    std::set<sgNodeID_t> fw_reached_nodes(sgNodeID_t n, int radius) const;

    /** @brief Returns all the nodes with connections in the graph
     *
     * @return set of node ids
     */
    std::unordered_set<sgNodeID_t>  get_connected_nodes() const;

    /** @brief Removes transitive connections form the graph
     *
     * @param radius
     */
    void remove_transitive_links(int radius);

    /**
     *
     */
    void report_connectivity();
    //void solve();

    bool are_connected(sgNodeID_t n1, sgNodeID_t n2);

    std::vector<std::vector<sgNodeID_t>> get_all_lines(uint16_t min_nodes, uint64_t min_total_size=0) const;
    std::vector<std::pair<sgNodeID_t,sgNodeID_t>> find_bubbles(uint32_t min_size,uint32_t max_size) ;

    void dump_to_text(std::string filename);
    void load_from_text(std::string filename);

    SequenceGraph & sg;
    std::vector<std::vector<Link>> links;

};
#endif //BSG_LINKAGEDIGRAPH_HPP
