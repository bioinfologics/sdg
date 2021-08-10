//
// Created by Bernardo Clavijo (EI) on 18/10/2017.
//

#ifndef SG_SEQUENCEGRAPH_HPP
#define SG_SEQUENCEGRAPH_HPP

#include <algorithm>
#include <array>
#include <vector>
#include <string>
#include <iostream>
#include <unordered_map>
#include <unordered_set>
#include <iosfwd>
#include <sdglib/indexers/UniqueKmerIndex.hpp>
#include "DistanceGraph.hpp"

class SequenceDistanceGraphPath;
class WorkSpace;
/**
 *
 * The most important features of this class are:
 *      - Nodes start from 1.
 *      - Nodes are signed to represent the direction in which they are being considered.
 */
class SequenceDistanceGraph : public DistanceGraph {
public:

    explicit SequenceDistanceGraph(WorkSpace & _ws):nodes(0),DistanceGraph(*this,"SDG"),ws(_ws) { //sdg gets initialised through LDG
        add_node(Node("",NodeStatus::Deleted)); //an empty deleted node on 0, just to skip the space
    };
    //SequenceDistanceGraph(const SequenceDistanceGraph &sg) = delete; // Avoid implicit generation of the copy constructor.

    bool operator==(const SequenceDistanceGraph &o) const {
        return (nodes == o.nodes and links==o.links);
    }

    friend std::ostream& operator<<(std::ostream &os, const SequenceDistanceGraph& sdg);


    std::string ls(int level=0,bool recursive=true) const;

    //=== I/O functions ===

    /**
     * @brief GFA loading function, it detects the format of the GFA file (1,2) and loads it appropriately
     * @param filename Path of the gfa file to load
     */
    void load_from_gfa(std::string filename);
    void load_from_bcalm(std::string filename,uint16_t k);

    void load_from_gfa1(std::ifstream &gfaf, std::ifstream &fastaf);
    void load_from_gfa2(std::ifstream &gfaf, std::ifstream &fastaf);

    void load_from_fasta(std::string filename);
    //TODO: move to DistanceGraph

    void write(std::ofstream & output_file);
    void read(std::ifstream & input_file);

    //=== read operations ===

    std::string get_node_sequence (sgNodeID_t n) const;
    uint64_t get_node_size (sgNodeID_t n) const;


    //=== graph operations ===
    /**
     * Adds a new node to the graph
     * @param n Node object to add
     * @return
     * Returns the ID of the added node
     */
    sgNodeID_t add_node(Node n);

    /**
     * Adds a new node to the graph from a string
     * @param n Node object to add
     * @return
     * Returns the ID of the added node
     */
    sgNodeID_t add_node(std::string seq, bool make_canonical=true);


    /**
     * Graph sanity check, makes sure the graph abides to the expected structure
     * @return
     * Whether the graph is valid or not
     */
    bool is_sane() const;

    /* TODO: deprecate and reimplement in DistanceGraph
     * Connected components, (TODO) optionally breaking up in repeats, nodes that class as repeats will be returned on their own
     */
    std::vector<std::vector<sgNodeID_t>> connected_components (int max_nr_totalinks=0, int max_nr_dirlinks=0, int min_rsize=0); //TODO: --> enable extra breaks in repeats


    /**
     * Delete a node
     * @param n ID of the node
     */
    void remove_node(sgNodeID_t n);
    //These two need to mark expanded edges, and transfer read maps and unique kmers for non-expanded, but just read map for expanded.

    /**
     * This creates a new node with the sequence of the full path, and connects to the same end connections as path.
     * Optionally removes the nodes that only participate in this path.
     * @param p
     * @param consume_nodes
     */
    void join_path(const SequenceDistanceGraphPath p, bool consume_nodes=true);
    // expand_path --> creates an edge with the consensus of a path, eliminates old nodes if only in path and unused edges

    uint32_t join_all_unitigs();
    //TODO: deprecate and replace/merge with get_all_lines
    std::vector<SequenceDistanceGraphPath> get_all_unitigs(uint16_t min_nodes);;

    /**
     * Makes multiple copies of a node to expand as repeat, connects to the bw and fw as specified,
     * removes those connections from the original node.
     *
     * @param nodeID
     * @param bw
     * @param fw
     */
    void expand_node(sgNodeID_t nodeID, std::vector<std::vector<sgNodeID_t>> bw,
                                    std::vector<std::vector<sgNodeID_t>> fw);

    /**
     * From a list of names get a list of node IDs from the graph
     * @param _oldnames String containing old names
     * @return
     * List of node IDs from the graph
     */
    std::vector<sgNodeID_t> oldnames_to_nodes(std::string _oldnames);

    //    std::vector<sgNodeID_t > find_canonical_repeats();

    /**
     * Get Name of a node
     * @param id Node ID in the graph
     * @return
     * Name of the node
     */
    const std::string& nodeID_to_name(sgNodeID_t id) const {
        return oldnames[std::abs(id)];
    }

    size_t count_active_nodes() const;

    void print_status();

    //=== internal variables ===

    std::vector<Node> nodes={};    /// Contains the actual nodes from the graph, nodes are generally accesed using its IDs on to this structure.
    std::string filename,fasta_filename;    /// Name of the files containing the graph and the fasta.
    std::vector<std::string> oldnames;      /// Mapping structure IDs to input names
    std::unordered_map<std::string,sgNodeID_t> oldnames_to_ids; /// Mapping structure from input names -> IDs

    WorkSpace &ws;

};
#endif //SG_SEQUENCEGRAPH_HPP
