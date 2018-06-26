//
// Created by Luis Yanes (EI) on 12/02/2018.
//

#ifndef BSG_LONGREADMAPPER_H
#define BSG_LONGREADMAPPER_H


#include <iostream>
#include <sglib/factories/KMerIDXFactory.h>
#include <sglib/readers/FileReader.h>
#include <sglib/readers/SequenceGraphReader.h>
#include <sglib/SMR.h>
#include <sglib/factories/StrandedMinSketchFactory.h>
#include <sglib/datastores/LongReadsDatastore.hpp>
#include <sglib/types/MappingTypes.hpp>
#include <minimap.h>

/**
 * Long read mapping to the graph, this class manages storage and computation of the alignments.
 */
class LongReadMapper {
    SequenceGraph & sg;
    mm_idx_t *graph_index = nullptr;
    mm_mapopt_t opt;
    uint8_t k=15;
    uint8_t w=10;

    /**
     * Stores an index of the mappings of a node to all the mappings where it appears.
     * This index can be queried to get information about all reads that map to a node.
     */
    std::vector< std::vector < std::vector<LongReadMapping>::size_type > > mappings_in_node;        /// Mappings matching node

    void update_indexes_from_mappings();
    LongReadMapping createMapping(uint32_t readID, const mm_reg1_t *regs0, int j, long long int node) const;
    void printMatch(const mm_idx_t *mi, std::ofstream &matchOutput, uint32_t readID, const std::string &read_name,
                    int read_len, const mm_reg1_t *regs0, int j);

public:

    LongReadMapper(SequenceGraph &sg, LongReadsDatastore &ds, uint8_t k=15, uint8_t w=10);
    ~LongReadMapper();

    void update_graph_index();
    LongReadsDatastore& getLongReadsDatastore() {return datastore;}


    void map_reads(std::unordered_set<uint32_t> readIDs = {});

    void read(std::string filename);

    void read(std::ifstream &inf);

    void write(std::string filename);

    void write(std::ofstream &ofs);



    LongReadsDatastore datastore;
    /**
     * This public member stores a flat list of mappings from the reads, it is accesed using the mappings_in_node index
     * or the read_to_mappings index.
     */
    std::vector<LongReadMapping> mappings;
    /**
     * Stores an index of the resulting mappings of a single long read, for each long read, stores the position of it's mappings.
     * This index can be used to query all the nodes that map to a single read.
     */
    std::vector< std::vector < std::vector<LongReadMapping>::size_type > > read_to_mappings;    /// Nodes in the read, 0 or empty = unmapped
};


#endif //BSG_LONGREADMAPPER_H
