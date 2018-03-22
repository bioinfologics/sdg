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
    std::vector<Link> get_bw_links( sgNodeID_t n);

    /*
     * Connected components, (TODO) optionally breaking up in repeats, nodes that class as repeats will be returned on their own
     */
    std::vector<std::vector<sgNodeID_t>> connected_components (int max_nr_totalinks=0, int max_nr_dirlinks=0, int min_rsize=0); //TODO: --> enable extra breaks in repeats

    // find bubbles in component of graph
    std::vector<std::vector<sgNodeID_t >> find_bubbles(std::vector<sgNodeID_t>);

    std::vector<nodeVisitor> depth_first_search(nodeVisitor node, unsigned int size_limit, unsigned int edge_limit, std::set<nodeVisitor> tabu={});

    std::vector<sgNodeID_t> breath_first_search(std::vector<sgNodeID_t> &nodes, unsigned int size_limit);

    SequenceGraphPath find_path_between(const sgNodeID_t seed, const sgNodeID_t target, unsigned int size_limit = 0, unsigned int edge_limit = 0);

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
    bool is_loop(std::array<sgNodeID_t, 4> nodes) {
        for (auto j = 0; j < 3; ++j)
            for (auto i = j + 1; i < 4; ++i)
                if (nodes[i] == nodes[j] or nodes[i] == -nodes[j]) return true; //looping node
        return false;
    }

    bool link_exists(sgNodeID_t from, sgNodeID_t to) {
        // Look for link between starting node and the new node.
        auto l = links[std::abs(from)].begin();
        // TODO: Can this just be a std::find?
        for (; l != links[std::abs(from)].end(); ++l){
            if (l->source == from and l->dest == to) break;
        }
        return l != links[std::abs(from)].end();
    }

    std::string& nodeID_to_name(sgNodeID_t id) {
        return oldnames[id];
    }

};

#endif //SG_SEQUENCEGRAPH_HPP
