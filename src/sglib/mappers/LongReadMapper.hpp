//
// Created by Luis Yanes (EI) on 12/02/2018.
//

#ifndef BSG_LONGREADMAPPER_H
#define BSG_LONGREADMAPPER_H


#include <iostream>
#include <sglib/factories/KMerIDXFactory.h>
#include <sglib/datastores/LongReadsDatastore.hpp>
#include <sglib/graph/SequenceGraph.hpp>
#include <sglib/types/MappingTypes.hpp>
#include <sglib/indexers/NKmerIndex.hpp>

/**
 * Long read mapping to the graph, this class manages storage and computation of the alignments.
 */
class LongReadMapper {

    const SequenceGraph & sg;

    uint8_t k=13;//15
    int window_size = 1500;//750;
    int window_slide = window_size/3;
    int min_score = 50;//11;
    int second_best_score_pct=90;
    int max_num_score_nodes = 100; //how many high-scoring nodes to consider per read

    /**
     * Stores an index of the mappings of a node to all the mappings where it appears.
     * This index can be queried to get information about all reads that map to a node.
     */
    std::vector< std::vector < std::vector<LongReadMapping>::size_type > > reads_in_node;        /// Reads matching node

    NKmerIndex assembly_kmers;

    void update_indexes_from_mappings();

public:

    LongReadMapper(SequenceGraph &sg, LongReadsDatastore &ds, uint8_t k=15);
    ~LongReadMapper();

    LongReadMapper operator=(const LongReadMapper &other);

    LongReadsDatastore& getLongReadsDatastore() {return datastore;}

    std::vector<uint64_t> get_node_read_ids(sgNodeID_t nodeID) const ;

    void map_reads(std::unordered_set<uint32_t> readIDs = {});

    void read(std::string filename);

    void read(std::ifstream &inf);

    void write(std::string filename);

    void write(std::ofstream &output_file);

    void update_graph_index();

    LongReadsDatastore datastore;
    /**
     * This public member stores a flat list of mappings from the reads, it is accessed using the mappings_in_node index
     * or the read_to_mappings index.
     */
    std::vector<LongReadMapping> mappings;
    /**
     * Stores an index of the resulting mappings of a single long read, for each long read, stores the position of it's mappings.
     * This index can be used to query all the nodes that map to a single read.
     */
    std::vector< std::vector < std::vector<LongReadMapping>::size_type > > read_to_mappings;    /// Nodes in the read, 0 or empty = unmapped

    static const bsgVersion_t min_compat;

};


#endif //BSG_LONGREADMAPPER_H
