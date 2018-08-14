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
#include <unordered_set>
#include <set>
#include "sglib/readers/FileReader.h"

typedef int64_t sgNodeID_t; //first node is 1; negatives are RC

#ifndef __hash128
#define __hash128
namespace std {
    //TODO: this hashing sucks, but it is needed
    template <> struct hash<__int128 unsigned>
    {
        size_t operator()(const __int128 unsigned & x) const
        {
            return hash<uint64_t>()((uint64_t)x);
        }
    };
}
#endif

enum sgNodeStatus_t {sgNodeActive,sgNodeDeleted};
struct graphPosition{
    sgNodeID_t node;
    uint32_t pos;
};
class Node{
public:
    Node(std::string _seq, sgNodeStatus_t _st) : sequence(_seq),status(_st){};
    Node(std::string _seq) : sequence(_seq),status(sgNodeActive){};
    bool operator==(const Node &o) const {
        return std::tie(status,sequence) == std::tie(o.status,o.sequence);
    }
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

class SequenceGraphPath;
class SequenceSubGraph;
class SequenceGraph {
public:
    SequenceGraph(){add_node(Node("",sgNodeStatus_t::sgNodeDeleted));};
    bool operator==(const SequenceGraph &o) const {
        return nodes == o.nodes;
    }
    bool is_sane();
    //=== I/O functions ===
    void load_from_gfa(std::string filename);
    void write_to_gfa(std::string filename,const std::unordered_set<sgNodeID_t> & marked_red={}, const std::vector<double> & depths={},
                      const std::unordered_set<sgNodeID_t> & selected_nodes={}, const std::vector<std::vector<Link>> & arg_links= {});
    void write(std::ofstream & output_file);
    void read(std::ifstream & input_file);

    //=== graph operations ===
    sgNodeID_t add_node(Node n);
    void add_link( sgNodeID_t source, sgNodeID_t dest, int32_t d);
    Link get_link( sgNodeID_t source, sgNodeID_t dest);
    std::vector<Link> get_fw_links( sgNodeID_t n);
    inline std::vector<Link> get_bw_links( sgNodeID_t n){ return get_fw_links (-n); };
    size_t count_active_nodes();

    /*
     * Connected components, (TODO) optionally breaking up in repeats, nodes that class as repeats will be returned on their own
     */
    std::vector<std::vector<sgNodeID_t>> connected_components (int max_nr_totalinks=0, int max_nr_dirlinks=0, int min_rsize=0); //TODO: --> enable extra breaks in repeats

    // find bubbles in component of graph
    std::vector<std::vector<sgNodeID_t >> find_bubbles(std::vector<sgNodeID_t>);
    std::vector<SequenceSubGraph> get_all_bubbly_subgraphs(uint32_t maxsubgraphs=0);
    void print_bubbly_subgraph_stats(const std::vector<SequenceSubGraph> & bubbly_paths);
    std::vector<std::pair<sgNodeID_t,int64_t>> get_distances_to(sgNodeID_t n, std::set<sgNodeID_t> destinations, int64_t max_dist);
    std::vector<SequenceGraphPath> find_all_paths_between(sgNodeID_t from,sgNodeID_t to, int64_t max_size, int max_nodes=20, bool abort_on_loops=true);

    // remove_node
    void remove_node(sgNodeID_t);
    // remove_link
    void remove_link(sgNodeID_t source, sgNodeID_t dest);
    //These two need to mark expanded edges, and transfer read maps and unique kmers for non-expanded, but just read map for expanded.

    /**
     * @brief expands a node creating as many copies as needed, then, distributes input and output links as per bw and fw
     * @param nodeID
     * @param bw
     * @param fw
     */
    void expand_node(sgNodeID_t nodeID, std::vector<std::vector<sgNodeID_t>> bw, std::vector<std::vector<sgNodeID_t>> fw);
    /**
     * @brief copies every middle node in the path, connects them, and moves connections from first and last
     * @param p
     */
    void expand_path(SequenceGraphPath p);
    void join_path(const SequenceGraphPath p, bool consume_nodes=true);
    // expand_path --> creates an edge with the consensus of a path, eliminates old nodes if only in path and unused edges
    uint32_t join_all_unitigs();
    std::vector<SequenceGraphPath> get_all_unitigs(uint16_t min_nodes);
    // simplify --> executes expand_path on every multi-sequence unitig
    std::vector<SequenceSubGraph> get_all_tribbles();
    void create_index(bool verbose=true);
    void create_63mer_index(bool verbose=true);
    void clear_index();

    // tip_clip -> eliminates tips.


    std::vector<sgNodeID_t> depth_first_search(const sgNodeID_t seed, unsigned int size_limit = 0, unsigned int edge_limit = 0, std::set<sgNodeID_t> tabu = {});

    //void explode_node( sgNodeID_t node, uint16_t k);
    //void explode_all_nodes ();
    //void collapse_identical_nodes ();

    //later
    //project spectra and use for flow


    //=== internal variables ===
    std::unordered_map<uint64_t, graphPosition> kmer_to_graphposition;
    std::unordered_map<__uint128_t, graphPosition> k63mer_to_graphposition;
    std::vector<Node> nodes;
    std::vector<std::vector<Link>> links;
    std::string filename,fasta_filename;
    std::vector<sgNodeID_t> oldnames_to_nodes(std::string _oldnames);
    std::vector<std::string> oldnames;
    std::unordered_map<std::string,sgNodeID_t> oldnames_to_ids;

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
    bool operator==(const SequenceGraphPath & other) const { return this->nodes==other.nodes; };
    SequenceGraphPath (const SequenceGraphPath& other): sg(other.sg),nodes(other.nodes){};
    SequenceGraphPath& operator=(const SequenceGraphPath& other){nodes=other.nodes;return *this;};
    explicit SequenceGraphPath(const SequenceGraph & _sg, const std::vector<sgNodeID_t> _nodes={})  : sg(_sg) ,nodes(_nodes) {};
    std::string get_fasta_header();
    std::string get_sequence();
    size_t get_sequence_size_fast();
    std::vector<Link> get_next_links() { return sg.get_fw_links(nodes.back());}
    bool extend_if_coherent(SequenceGraphPath s);
    void reverse();
    bool is_canonical();
    bool is_unitig();
    const bool operator< (const SequenceGraphPath & other);
    const bool operator== (const SequenceGraphPath & other);

private:
    const SequenceGraph& sg;
};


class SequenceSubGraph {
public:
    std::vector<sgNodeID_t> nodes;
    explicit SequenceSubGraph(SequenceGraph & _sg, std::vector<sgNodeID_t> _nodes={})  : sg(_sg) ,nodes(_nodes) {};
    SequenceGraphPath make_path(); //returns empty path if not linear
    const uint64_t total_size () const;
private:
    SequenceGraph& sg;
};

#endif //SG_SEQUENCEGRAPH_HPP
