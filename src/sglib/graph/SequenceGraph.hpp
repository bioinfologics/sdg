//
// Created by Bernardo Clavijo (EI) on 18/10/2017.
//

#ifndef SG_SEQUENCEGRAPH_HPP
#define SG_SEQUENCEGRAPH_HPP

#include <vector>
#include <algorithm>
#include <string>
#include <map>
#include <unordered_map>
#include <set>
#include <iostream>
#include <array>
#include <unordered_set>
#include <sglib/types/KmerTypes.hpp>
#include <sglib/types/GenericTypes.hpp>
#include <sglib/graph/SequenceSubGraph.hpp>
#include <sglib/graph/SequenceGraphPath.hpp>
#include <sglib/logger/OutputLog.h>

class SequenceGraphPath;
class SequenceSubGraph;
/**
 * @brief Class representing sequence graphs
 *
 * This class represents the assembly problem as a flow formulation.
 * The most important features of this class are:
 *      - Nodes start from 1.
 *      - Nodes are signed to represent the direction in which they are being considered.
 *      - Nodes are saved in canonical form.
 */
class SequenceGraph {
public:
    //=== internal variables ===

    std::vector<Node> nodes;    /// Contains the actual nodes from the graph, nodes are generally accesed using its IDs on to this structure.
    std::vector<std::vector<Link>> links;   /// List of all links, links are stored in their canonical from (-1 -> 2)
    std::string filename,fasta_filename;    /// Name of the files containing the graph and the fasta.
    std::vector<std::string> oldnames;      /// Mapping structure IDs to input names
    std::unordered_map<std::string,sgNodeID_t> oldnames_to_ids; /// Mapping structure from input names -> IDs
    std::unordered_map<uint64_t, graphPosition> kmer_to_graphposition;  /// Indexing structure saves a unique kmer -> graphPosition
    std::unordered_map<__uint128_t, graphPosition> k63mer_to_graphposition; /// Indexing structure saves a unique 63mer -> graphPosition

    bool operator==(const SequenceGraph &o) const {
        return nodes == o.nodes;
    }

    SequenceGraph(){};
    SequenceGraph(const SequenceGraph &sg) = delete; // Avoid implicit generation of the copy constructor.
    //=== I/O functions ===
    void load_from_gfa(std::string filename);
    void write_to_gfa(std::string filename, const std::vector<std::vector<Link>> &arg_links={},
                          const std::vector<sgNodeID_t> &selected_nodes={}, const std::unordered_set<sgNodeID_t> &mark_red={},
                          const std::vector<double> &depths={});
    void write(std::ofstream & output_file);
    void read(std::ifstream & input_file);

    //=== graph operations ===
    /**
     * Adds a new node to the graph
     * @param n Node object to add
     * @return
     * Returns the ID of the added node
     */
    sgNodeID_t add_node(Node n);
    void add_link( sgNodeID_t source, sgNodeID_t dest, int32_t d);

    /**
     *
     * @param source
     * @param dest
     * @return
     */
    Link get_link(sgNodeID_t source, sgNodeID_t dest);

    /**
     * Get list of Links from a node
     * @param n Node to search links for
     * @return
     * Returns a vector of Links leaving from n
     */
    std::vector<Link> get_fw_links( sgNodeID_t n) const ;
    /**
     * Get list of links to a node
     * @param n Node to search links for
     * @return
     * Returns a vector of Links going to n
     */
    std::vector<Link> get_bw_links( sgNodeID_t n) const ;

    /**
     * Find node IDs for forward nodes of n
     * @param n Node to search neighbours for
     * @return
     * Returns a vector of forward node IDs
     */
    std::vector<sgNodeID_t> get_fw_nodes(sgNodeID_t n) const;
    /**
     * Find node IDs for backward nodes of n
     * @param n Node to search neighbours for
     * @return
     * Returns a vector of backward node IDs
     */
    std::vector<sgNodeID_t> get_bw_nodes(sgNodeID_t n) const;

    /**
     * Find node IDs for neighbouring nodes of n
     * @param n Node to search neighbours for
     * @return
     * Returns a vector of neighbouring node IDs
     */
    std::vector<sgNodeID_t> get_neighbour_nodes(sgNodeID_t n) const {
        std::unordered_set<sgNodeID_t > result;
        auto fwns(get_fw_nodes(n));
        auto bwns(get_bw_nodes(n));
        for (const auto &fn:fwns) {
            result.insert(fn);
        }
        for (const auto &bn:bwns) {
            result.insert(bn);
        }
        return std::vector<sgNodeID_t > (result.begin(), result.end());
    }

    /**
     * Graph sanity check, makes sure the graph abides to the expected structure
     * @return
     * Whether the graph is valid or not
     */
    bool is_sane() const;

    /*
     * Connected components, (TODO) optionally breaking up in repeats, nodes that class as repeats will be returned on their own
     */
    std::vector<std::vector<sgNodeID_t>> connected_components (int max_nr_totalinks=0, int max_nr_dirlinks=0, int min_rsize=0); //TODO: --> enable extra breaks in repeats

    // find bubbles in component of graph
    std::vector<std::vector<sgNodeID_t >> find_bubbles(std::vector<sgNodeID_t>);

    /**
     * This function returns all nodes that can be reached again if following the links of the next
     * _complexity_ nodes, in particular, this function returns the nodes represented
     * by ****** in the following image:
     *
     *\verbatim
     *            ^^^^^^^^^^^
     *           |           |
     *  ========---*********---========
     *\endverbatim
     *  The number of elements in the outer (^^^^^) part of the loop is defined by "complexity"
     * @param complexity Maximum number of nodes in the loop complexity
     * @return IDs of all the nodes involved in loops
     *
    */
    std::vector<sgNodeID_t > get_loopy_nodes(int complexity=3);

    /**
     *
     * @param loopy_node
     * @return IDs of flanking nodes to a repeat
     * Returns the ========= nodes from get_loopy_nodes when passed the ********** node
     */
    std::vector<sgNodeID_t> get_flanking_nodes(sgNodeID_t loopy_node);

    /**
     * Finds all reachable nodes within size_limit and edge_limit from the specified node, returns a vector with edge and base distance
     * foreach reachable node.
     * @param node Node to start exploring
     * @param size_limit Range of bases within search, i.e if distance between node
     * @param edge_limit Number of nodes in a single path to consider
     * @param tabu List of nodes not to visit
     * @return Returns a vector of nodeVisitors storing the path distance and bases distance between origin of search and each node
     */
    std::vector<nodeVisitor> depth_first_search (nodeVisitor node, unsigned int size_limit, unsigned int edge_limit, std::set<nodeVisitor> tabu={});

    std::vector<nodeVisitor> breath_first_search(nodeVisitor seed, unsigned int size_limit, unsigned int edge_limit, std::set<nodeVisitor> tabu={});

    /**
     * Finds a path between two nodes within size_limit and edge_limit.
     * If any of the limits is violated before finding a path, an empty path is returned
     * @param seed
     * @param target
     * @param size_limit
     * @param edge_limit
     * @return
     */
    SequenceGraphPath find_path_between(const sgNodeID_t seed, const sgNodeID_t target, unsigned int size_limit = 0, unsigned int edge_limit = 0);

    /**
     * Explores around "nodes" in both directions until either the size_limit or the edge_limit is reached,
     * considers all nodes multiple times always taking the shortest path as the limit (I.E. if there's a repeat,
     * the shortest path to the repeated node is taken).
     * @param nodes Nodes to start exploring from
     * @param size_limit Limit in number of bases to explore
     * @param edge_limit Limit in number of edges to go through (length of the path until a node)
     * @return Returns a vector of all nodes in the solution
     *
     */
    std::vector<nodeVisitor> explore_nodes(std::vector<std::string> &nodes, unsigned int size_limit, unsigned int edge_limit);

    /**
     * Delete a node
     * @param n ID of the node
     */
    void remove_node(sgNodeID_t n);
    // remove_link
    /**
     * Delete a link
     * @param source source
     * @param dest  destination
     * Returns true if a link was removed, false otherwise
     */
    bool remove_link(sgNodeID_t source, sgNodeID_t dest);
    //These two need to mark expanded edges, and transfer read maps and unique kmers for non-expanded, but just read map for expanded.

    void join_path(const SequenceGraphPath p, bool consume_nodes=true);
    // expand_path --> creates an edge with the consensus of a path, eliminates old nodes if only in path and unused edges
    uint32_t join_all_unitigs();
    std::vector<SequenceGraphPath> get_all_unitigs(uint16_t min_nodes);
    std::vector<SequenceSubGraph> get_all_tribbles(){};


    void expand_node(sgNodeID_t nodeID, std::vector<std::vector<sgNodeID_t>> bw,
                                    std::vector<std::vector<sgNodeID_t>> fw);
    std::vector<std::pair<sgNodeID_t,int64_t>> get_distances_to(sgNodeID_t n, std::set<sgNodeID_t> destinations, int64_t max_dist);
    std::vector<SequenceSubGraph> get_all_bubbly_subgraphs(uint32_t maxsubgraphs = 0);
    void print_bubbly_subgraph_stats(const std::vector<SequenceSubGraph> &bubbly_paths);

    /**
     * From a list of names get a list of node IDs from the graph
     * @param _oldnames String containing old names
     * @return
     * List of node IDs from the graph
     */
    std::vector<sgNodeID_t> oldnames_to_nodes(std::string _oldnames);

    /**
     * Get the node object containing all node information
     * @param n Node ID from the graph
     * @return
     * Node sequence and status
     */
    Node& get_node(sgNodeID_t n) { return nodes[(n>0)?n:-n];}

//    std::vector<sgNodeID_t > find_canonical_repeats();

    bool is_loop(std::array<sgNodeID_t, 4> nodes);

    /**
     * Checks existance of a link in the graph
     * @param from  Directed node
     * @param to    Directed node
     * @return
     * True if the Link on the direction of the nodes provided exists
     */
    bool link_exists(sgNodeID_t from, sgNodeID_t to) const {
        // Look for link between starting node and the new node.
        auto l = links[std::abs(from)].begin();
        // TODO: Can this just be a std::find?
        for (; l != links[std::abs(from)].end(); ++l){
            if (l->source == from and l->dest == to) break;
        }
        return l != links[std::abs(from)].end();
    }

    /**
     * Get Name of a node
     * @param id Node ID in the graph
     * @return
     * Name of the node
     */
    const std::string& nodeID_to_name(sgNodeID_t id) const {
        return oldnames[std::abs(id)];
    }

    /**
     * Function to index the graph
     * Stores the result in the local kmer_to_graphposition object
     */
    void create_index(bool verbose=true);

    /**
     * Function to index the graph
     * Stores the result in the local kmer_to_graphposition object
     */
    void create_63mer_index(bool verbose=true);


    size_t count_active_nodes();
    std::vector<SequenceGraphPath> find_all_paths_between(sgNodeID_t from,sgNodeID_t to, int64_t max_size, int max_nodes=20, bool abort_on_loops=true);};

#endif //SG_SEQUENCEGRAPH_HPP
