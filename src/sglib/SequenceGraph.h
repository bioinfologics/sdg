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

typedef int64_t sgNodeID_t; //first node is 1; negatives are RC

enum sgNodeStatus_t {sgNodeActive,sgNodeDeleted};

class Node{
public:
    Node(std::string _seq) : sequence(_seq),status(sgNodeActive){};
    std::string sequence;
    sgNodeStatus_t status;
    bool is_canonical();
    void make_rc();
};

class Link{
public:
    Link( sgNodeID_t _src, sgNodeID_t _dst, int32_t _dist) : source(_src), dest(_dst), dist(_dist) {};
    sgNodeID_t source,dest;
    int32_t dist;

    bool operator==( const  Link);
    bool operator<(const Link)const;

    friend std::ostream &operator<<(std::ostream &os, const Link &link) {
        os << link.source << " -> " << link.dest;
        return os;
    }

};

struct nodeVisitor {
    sgNodeID_t node;
    uint dist;
    uint path_length;
    nodeVisitor(sgNodeID_t n, uint d, uint p) : node(n), dist(d), path_length(p) {}
    bool operator<(const nodeVisitor &o) const {return node < o.node;}
    bool operator==(const nodeVisitor &o) const {return node == o.node;}
    nodeVisitor reverseDirection() const {
        return nodeVisitor(-node, dist, path_length);
    }
};

class SequenceGraphPath;
class SequenceSubGraph;
class SequenceGraph {
public:
    SequenceGraph(){};
    //=== I/O functions ===
    void load_from_gfa(std::string filename);
    void write_to_gfa(std::string filename);

    //=== graph operations ===
    sgNodeID_t add_node(Node n);
    void add_link( sgNodeID_t source, sgNodeID_t dest, int32_t d);

    std::vector<Link> get_fw_links( sgNodeID_t n);
    std::vector<Link> get_bw_links( sgNodeID_t n);

    /*
     * Connected components, (TODO) optionally breaking up in repeats, nodes that class as repeats will be returned on their own
     */
    std::vector<std::vector<sgNodeID_t>> connected_components (int max_nr_totalinks=0, int max_nr_dirlinks=0, int min_rsize=0); //TODO: --> enable extra breaks in repeats

    // find bubbles in component of graph
    std::vector<std::vector<sgNodeID_t >> find_bubbles(std::vector<sgNodeID_t>);

    std::vector<nodeVisitor> depth_first_search(nodeVisitor node, unsigned int size_limit, unsigned int edge_limit, std::set<nodeVisitor> tabu={});

    std::vector<sgNodeID_t> breath_first_search(std::vector<sgNodeID_t> &nodes, unsigned int size_limit);

    // remove_node
    void remove_node(sgNodeID_t);
    // remove_link
    void remove_link(sgNodeID_t source, sgNodeID_t dest);
    //These two need to mark expanded edges, and transfer read maps and unique kmers for non-expanded, but just read map for expanded.

    void join_path(const SequenceGraphPath p, bool consume_nodes=true);
    // expand_path --> creates an edge with the consensus of a path, eliminates old nodes if only in path and unused edges
    void join_all_unitigs();
    std::vector<SequenceGraphPath> get_all_unitigs(uint16_t min_nodes);
    // simplify --> executes expand_path on every multi-sequence unitig


    // tip_clip -> eliminates tips.


    //void explode_node( sgNodeID_t node, uint16_t k);
    //void explode_all_nodes ();
    //void collapse_identical_nodes ();

    //later
    //project spectra and use for flow


    //=== internal variables ===
    std::vector<sgNodeID_t> oldnames_to_nodes(std::string _oldnames);
    std::vector<Node> nodes;
    std::vector<std::vector<Link>> links;
    std::unordered_map<std::string,sgNodeID_t> oldnames_to_ids;
    std::vector<std::string> oldnames;
    std::string filename,fasta_filename;

    void consume_nodes(const SequenceGraphPath &p, const std::set<sgNodeID_t> &pnodes);

    std::vector<sgNodeID_t > find_canonical_repeats();
    bool is_loop(std::array<sgNodeID_t, 4> nodes) {
        for (auto j = 0; j < 3; ++j)
            for (auto i = j + 1; i < 4; ++i)
                if (nodes[i] == nodes[j] or nodes[i] == -nodes[j]) return true; //looping node
        return false;
    }
};


class SequenceGraphPath {
public:
    std::vector<sgNodeID_t> nodes;
    explicit SequenceGraphPath(SequenceGraph & _sg, const std::vector<sgNodeID_t> _nodes={})  : sg(_sg) ,nodes(_nodes) {};
    std::string get_fasta_header();
    std::string get_sequence();
    void reverse();
    bool is_canonical();

private:
    SequenceGraph& sg;
};

class SequenceSubGraph {
public:
    std::vector<sgNodeID_t> nodes;
    explicit SequenceSubGraph(SequenceGraph & _sg, std::vector<sgNodeID_t> _nodes={})  : sg(_sg) ,nodes(_nodes) {};
    SequenceGraphPath make_path(); //returns empty path if not linear

    void write_to_gfa(std::string filename);
private:
    SequenceGraph& sg;
};

#endif //SG_SEQUENCEGRAPH_HPP
