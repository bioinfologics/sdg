//
// Created by Bernardo Clavijo (EI) on 18/10/2017.
//

#ifndef SG_SEQUENCEGRAPH_HPP
#define SG_SEQUENCEGRAPH_HPP

#include <vector>
#include <string>
#include <unordered_map>

typedef int64_t sgNodeID_t; //first node is 1; negatives are RC



class Node{
public:
    Node(std::string _seq) : sequence(_seq){};
    std::string sequence;
    uint8_t status;
    bool is_canonical();
    void make_rc();
};

class Link{
public:
    Link( sgNodeID_t _src, sgNodeID_t _dst, int32_t _dist) : source(_src), dest(_dst), dist(_dist) {};
    sgNodeID_t source,dest;
    int32_t dist;
};

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
    inline std::vector<Link> get_bw_links( sgNodeID_t n){ return get_fw_links (-n); };

    /*
     * Connected components, optionally breaking up in repeats, nodes that class as repeats will be returned on their own
     */
    std::vector<std::vector<sgNodeID_t>> connected_components (int max_nr_totalinks=0, int max_nr_dirlinks=0, int min_rsize=0); //--> enable extra breaks in repeats

    // remove_node
    // remove_link

    //These two need to mark expanded edges, and transfer read maps and unique kmers for non-expanded, but just read map for expanded.
    // expand_path --> creates an edge with the consensus of a path, eliminates old nodes if only in path and unused edges
    // simplify --> executes expand_path on every multi-sequence unitig

    //find balanced repeat subgraphs


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

};


class SequenceGraphPath {
public:
    std::vector<sgNodeID_t> nodes;
    explicit SequenceGraphPath(SequenceGraph & _sg, std::vector<sgNodeID_t> _nodes={})  : sg(_sg) ,nodes(_nodes) {};
    std::string get_fasta_header();
    std::string get_sequence();

private:
    SequenceGraph& sg;
};




#endif //SG_SEQUENCEGRAPH_HPP
