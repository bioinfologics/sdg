//
// Created by Bernardo Clavijo (EI) on 12/05/2018.
//

#pragma once

#include <map>
#include <unordered_set>
#include <fstream>

#include "sdglib/types/MappingTypes.hpp"
#include "sdglib/factories/KMerIDXFactory.hpp"
#include <sdglib/indexers/UniqueKmerIndex.hpp>
#include <sdglib/Version.hpp>

class PairedReadsDatastore;
class UniqueKmerIndex;
class Unique63merIndex;
class WorkSpace;

/**
 * A mapper for linked reads from a PairedReadsDatastore.
 * Supports partial remapping of unmapped reads or of a selection list.
 */
class PairedReadsMapper {

public:
    PairedReadsMapper(WorkSpace &_ws, PairedReadsDatastore &_datastore);

    /**
     * @brief Provides an overview of the information in the PairedReadsMapper
     * @param level Base indentation level to use on the result
     * @param recursive Whether it should explore or not the rest of the hierarchy
     * @return
     * A text summary of the information contained in a PairedReadsMapper
     */
    std::string ls(int level=0,bool recursive=true) const;

    friend std::ostream& operator<<(std::ostream &os, const PairedReadsMapper &prm);

    void write(std::ofstream & output_file);
    void read(std::ifstream & input_file);
    void dump_readpaths(std::ofstream &opf);
    void load_readpaths(std::ifstream &ipf);

    /** @brief creates a read path for each read through mapping
     *
     * @return
     */
    void path_reads(uint8_t k=63,int filter=200, bool fill_offsets=false);

    /** @brief returns a list of possible nodes fw, node needs to appear with same sign in read's path, negative read id means reversed path
     *
     * @param read_id
     * @param node
     * @param use_pair
     * @param collapse_pair
     */
    std::vector<sgNodeID_t> path_fw(seqID_t read_id,sgNodeID_t node, bool use_pair=true, bool collapse_pair=true) const;

    std::vector<std::vector<sgNodeID_t>> all_paths_fw(sgNodeID_t node, bool use_pair=true, bool collapse_pair=true) const;

    /** @brief Estimates the fragment size distribution for the mapped data-store using the mapped pairs.
     *
     * Estimates the fragment size distribution for the mapped data-store using the mapped pairs.
     * The result is also stored in this.rfdist and this.frdist. and one of the vectors is also returned depending
     * on majority of the orientations of the reads
     *
     * @return distance count vector, the same as this.rfdist or this.frdist depending on the majority of the orientations of the reads
     */
    std::vector<uint64_t> size_distribution();

    /** @brief Populates the read_direction_in_node collection for the dataset
     *
     */
    void populate_orientation();
    /** @brief Prints the count of pairs mapped.
     *
     * Prints the count of pairs mapped where no end mapped, a single end mapped and both ends mapped and of
     * those how many mapped to a single node.
     */
    void print_status() const ;

    PairedReadsMapper& operator=(const PairedReadsMapper &other);

    /** Get the most conected neighbour of a node
     *
     */
    sgNodeID_t get_node_inmediate_neighbours(sgNodeID_t node);

    std::vector<int64_t> get_paths_in_node(sgNodeID_t nid);

    WorkSpace &ws;
    const PairedReadsDatastore & datastore;
    /**
     * reads_in_node[i] contains all the read ids from the data-store mapped to node i
     */
    std::vector<std::vector<ReadMapping>> reads_in_node;

    /**
     * read_to_node[i] contains the node where read i was mapped
     */
    std::vector<sgNodeID_t> read_to_node; //id of the main node if mapped, set to 0 to remap on next process

    /**
     * read_direction_in_node[i] has the direction with which read i was mapped in the corresponding mapping (in read_to_node[i])
     */
    std::vector<bool> read_direction_in_node;//0-> fw, 1->rev;

    /**
     * Reverse forward distance accumulator, do not use directly, if you need the fragment size use size_distribution()
     */
    std::vector<uint64_t> rfdist; /// Revese forward distance accumulator
    /**
     * Forward reverse distance accumulator, do not use directly, if you need the fragment size use size_distribution()
     */
    std::vector<uint64_t> frdist; /// Forward reverse distance accumulator

    std::vector<ReadPath> read_paths;

    std::vector<std::vector<std::pair<uint32_t,uint32_t>>> read_path_offsets;

    std::vector<std::vector<int64_t>> paths_in_node;
    static const sdgVersion_t min_compat;

};

