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
#include <sglib/readers/FileReader.h>
#include <sglib/graph/SequenceSubGraph.hpp>
#include <sglib/graph/SequenceGraphPath.hpp>

class SequenceGraphPath;
class SequenceSubGraph;
class SequenceGraph {
public:
    SequenceGraph(){};
    //=== I/O functions ===
    void load_from_gfa(std::string filename);
    void write_to_gfa(std::string filename,const std::unordered_set<sgNodeID_t> & marked_red={}, const std::vector<double> & depths={}, const std::unordered_set<sgNodeID_t> & selected_nodes={});
    void write(std::ofstream & output_file);
    void read(std::ifstream & input_file);

    //=== graph operations ===
    sgNodeID_t add_node(Node n);
    void add_link( sgNodeID_t source, sgNodeID_t dest, int32_t d);

    std::vector<Link> get_fw_links( sgNodeID_t n) const ;
    std::vector<Link> get_bw_links( sgNodeID_t n) const ;
    bool is_sane() const;

    /*
     * Connected components, (TODO) optionally breaking up in repeats, nodes that class as repeats will be returned on their own
     */
    std::vector<std::vector<sgNodeID_t>> connected_components (int max_nr_totalinks=0, int max_nr_dirlinks=0, int min_rsize=0); //TODO: --> enable extra breaks in repeats

    // find bubbles in component of graph
    std::vector<std::vector<sgNodeID_t >> find_bubbles(std::vector<sgNodeID_t>);

    /**
    * This function returns all nodes that can be reached again if following the links of the next
    * _complexity_ nodes
    * @param complexity Maximum number of nodes in the loop complexity
    * @return IDs of all the nodes involved in loops
    */
    std::vector<sgNodeID_t > get_loopy_nodes(uint complexity=3);
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
    std::vector<nodeVisitor> explore_nodes(std::vector<std::string> &nodes, uint size_limit, uint edge_limit);

    // remove_node
    void remove_node(sgNodeID_t);
    // remove_link
    void remove_link(sgNodeID_t source, sgNodeID_t dest);
    //These two need to mark expanded edges, and transfer read maps and unique kmers for non-expanded, but just read map for expanded.

    void join_path(const SequenceGraphPath p, bool consume_nodes=true);
    // expand_path --> creates an edge with the consensus of a path, eliminates old nodes if only in path and unused edges
    void join_all_unitigs();
    std::vector<SequenceGraphPath> get_all_unitigs(uint16_t min_nodes);
    std::vector<SequenceSubGraph> get_all_tribbles();


    void expand_node(sgNodeID_t nodeID, std::vector<std::vector<sgNodeID_t>> bw,
                                    std::vector<std::vector<sgNodeID_t>> fw);
    std::vector<std::pair<sgNodeID_t,int64_t>> get_distances_to(sgNodeID_t n, std::set<sgNodeID_t> destinations, int64_t max_dist);
    std::vector<SequenceSubGraph> get_all_bubbly_subgraphs(uint32_t maxsubgraphs = 0);
    std::vector<SequenceGraphPath> find_all_paths_between(sgNodeID_t from,sgNodeID_t to, int64_t max_size);
    // simplify --> executes expand_path on every multi-sequence unitig


    // tip_clip -> eliminates tips.


    //void explode_node( sgNodeID_t node, uint16_t k);
    //void explode_all_nodes ();
    //void collapse_identical_nodes ();

    //later
    //project spectra and use for flow


    //=== internal variables ===

    std::vector<Node> nodes;
    std::vector<std::vector<Link>> links;
    std::string filename,fasta_filename;
    std::vector<sgNodeID_t> oldnames_to_nodes(std::string _oldnames);
    std::vector<std::string> oldnames;
    std::unordered_map<std::string,sgNodeID_t> oldnames_to_ids;

    void consume_nodes(const SequenceGraphPath &p, const std::set<sgNodeID_t> &pnodes);

    std::vector<sgNodeID_t > find_canonical_repeats();

    bool is_loop(std::array<sgNodeID_t, 4> nodes);

    bool link_exists(sgNodeID_t from, sgNodeID_t to) const {
        // Look for link between starting node and the new node.
        auto l = links[std::abs(from)].begin();
        // TODO: Can this just be a std::find?
        for (; l != links[std::abs(from)].end(); ++l){
            if (l->source == from and l->dest == to) break;
        }
        return l != links[std::abs(from)].end();
    }

    const std::string& nodeID_to_name(sgNodeID_t id) const {
        return oldnames[id];
    }

};

#endif //SG_SEQUENCEGRAPH_HPP
