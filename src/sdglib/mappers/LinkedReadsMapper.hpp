//
// Created by Bernardo Clavijo (EI) on 11/02/2018.
//

#pragma once

#include <set>
#include <unordered_set>
#include <map>
#include <sdglib/types/MappingTypes.hpp>
#include <sdglib/indexers/UniqueKmerIndex.hpp>
#include <sdglib/Version.hpp>

class WorkSpace;
class UniqueKmerIndex;
class Unique63merIndex;
class LinkedReadsDatastore;

using LinkedTag = uint32_t;

class TagNeighbour {
public:
    TagNeighbour(){};
    TagNeighbour(sgNodeID_t n, float s):node(n),score(s){};
    sgNodeID_t node;
    float score; //breaking the latte principle
};

/**
 * @brief A mapper for linked reads from a LinkedReadsDatastore.
 *
 * Supports partial remapping of unmapped reads or of a selection list.
 */
class LinkedReadsMapper {

public:
    LinkedReadsMapper(WorkSpace &_ws, LinkedReadsDatastore &_datastore);

    /**
     * @brief Provides an overview of the information in the LinkedReadsMapper
     * @param level Base indentation level to use on the result
     * @param recursive Whether it should explore or not the rest of the hierarchy
     * @return
     * A text summary of the information contained in a LinkedReadsMapper
     */
    std::string ls(int level=0,bool recursive=true) const;

    friend std::ostream& operator<<(std::ostream &os, const LinkedReadsMapper &lirm);

    void write(std::ofstream & output_file);
    void read(std::ifstream & input_file);
    void dump_readpaths(std::ofstream &opf);
    void load_readpaths(std::ifstream &ipf);



    /**
     * Checks that the graph and datastore are the same and if they are assigns reads_in_node and nodes_in_read collections
     *
     * @param other LinkedReadsMapper to compare with
     * @return
     */
    LinkedReadsMapper& operator=(const LinkedReadsMapper &other);

    /** @brief creates a read path for each read through mapping
     *
     * @return
     */
    void path_reads(uint8_t k=63,int filter=200, bool fill_offsets=false);

    /** @brief Prints the count of pairs mapped.
     *
     * Prints the count of pairs mapped where no end mapped, a single end mapped and both ends mapped and of those how many mapped to a single node.
     */
    void print_status() const;

    /** @brief Given a nodeID returns a set of all tags mapped to that node.
     * If there are no mapped tags returns an empty set
     *
     * @param n nodeID
     * @return set of bsg10xtags
     */
    std::set<LinkedTag> get_node_tags(sgNodeID_t n);

    /**
     * Creates a nieghbours matrix with all nodes where for each node the function finds all nodes that have tags that
     * cover min_score of the total tags of each node.
     *
     * Example:
     *      - node A has reads with tags 5, 5, 6, 7, 8
     *      - node B has reads with tags 5, 8, 9, 9, 10
     *
     * Then B tags cover 3/6=0.5 of A reads, and A tags cover 2/5=0.4 of B reads.
     * If min_score=.5 then B is in A's neighbours, but A is not in B's
     *
     * Results are stored in the tag_neighbours vector
     * @param min_size
     * @param min_score
     * @param min_mapped_reads_per_tag
     */
    void compute_all_tag_neighbours(int min_size,float min_score, int min_mapped_reads_per_tag=2);

    void write_tag_neighbours(std::string filename);
    void read_tag_neighbours(std::string filename);

    WorkSpace &ws;
    LinkedReadsDatastore &datastore;

    /**
     * Collection of read mappings
     * reads_in_node[nodeID] = [vector of mappings to nodeID... ]
     */
    std::vector<std::vector<ReadMapping>> reads_in_node;

    /**
     * Read to node index
     * read_to_node[readID] = nodeID where the read is mapped to
     */
    std::vector<sgNodeID_t> read_to_node;//id of the main node if mapped, set to 0 to remap on next process

    /**
     *  Collection of neighbouring nodes
     *  Completed using compute_all_tag_neighbours()
     *  tag_neighbours[nodeID] = [collection of 10 determined neighbours (see TagNeighbour) ...]
     */
    std::vector<std::vector<TagNeighbour>> tag_neighbours; //not persisted yet!
    std::vector<ReadPath> read_paths;

    std::vector<std::vector<std::pair<uint32_t,uint32_t>>> read_path_offsets;

    std::vector<std::vector<int64_t>> paths_in_node;

    static const sdgVersion_t min_compat;
};
