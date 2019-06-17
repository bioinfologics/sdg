//
// Created by Bernardo Clavijo (EI) on 11/02/2018.
//

#pragma once

#include <map>

#include "sdglib/graph/SequenceDistanceGraph.hpp"
#include <sdglib/types/MappingTypes.hpp>
#include <sdglib/indexers/UniqueKmerIndex.hpp>
#include <sdglib/Version.hpp>
class WorkSpace;
class UniqueKmerIndex;
class Unique63merIndex;
class LinkedReadsDatastore;

using bsg10xTag = uint32_t;

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
    LinkedReadsMapper(const WorkSpace &_ws, LinkedReadsDatastore &_datastore);
    void write(std::ofstream & output_file);
    void read(std::ifstream & input_file);

    /** @brief Maps the provided LinkedReadsDatastore to the provided graph, saves the results to the reads_in_node and nodes_in_read collections.
     * The mapper only considers unique perfect matches, if a single kmer from the read matches 2 different nodes the read is marked as multi-mapping and the mapping is discarded
     * If a reads_to_remap set is passed to the function only the selected set of reads is mapped, the rest of the mappings remains as they are
     * TODO: readIDs to alias ?
     *
     * To access the results see the reads_in_node and nodes_in_read collections in this class
     *
     * @param reads_to_remap set of readsIDs to remap
     */
    void map_reads(std::unordered_set<uint64_t> const &  reads_to_remap={});

    /** @brief Clears the reads_in_node and nodes_in_read collections and runs the map_reads() function with no selected set (maps all reads again)
     *
     */
    void remap_all_reads();

    /**
     * Cheks that the graph and datastore are the same and if they are assigns reads_in_node and nodes_in_read collections
     *
     * @param other
     * @return
     */
    LinkedReadsMapper operator=(const LinkedReadsMapper &other);

    /**
     * This is the same as map_reads() but using k63
     * TODO: merge with map_reads()
     * @param reads_to_remap
     */
    void map_reads63(std::unordered_set<uint64_t> const &  reads_to_remap={});

    /**
     * Same as remap_all_reads but with k63
     * TODO: merge with remap_all_reads()
     */
    void remap_all_reads63();

    /** @brief Clears the mappings for the nodes marked as deleted (sg.nodes[n].status==sgNodeDeleted) in the graph
     *
     */
    void remove_obsolete_mappings();

    /*void remap_reads();
    uint64_t process_reads_from_file(uint8_t, uint16_t, std::unordered_map<uint64_t , graphPosition> &, std::string , uint64_t, bool tags=false, std::unordered_set<uint64_t> const & reads_to_remap={});
    void save_to_disk(std::string filename);
    void load_from_disk(std::string filename);*/

    /** @brief Prints the count of pairs mapped.
     *
     * Prints the count of pairs mapped where no end mapped, a single end mapped and both ends mapped and of
     * those how many mapped to a single node.
     */
    void print_status();

    /** @brief Given a nodeID returns a set of all tags mapped to that node.
     * if there are no mapped tags returns an empty set
     *
     * @param n nodeID
     * @return set of bsg10xtags
     */
    std::set<bsg10xTag> get_node_tags(sgNodeID_t n);

    /**
     * Returns a tags_to_nodes type map, a collection of tags with an associated vector of nodes where each tag mapped
     * Each tag can map to more than one node ( multiple reads mapping to different nodes with the same tag ) map[tag] = [node1, node2, noden...]
     * TODO: tags should be plural
     * @param min_nodes minimum ammount of nodes in the vector for a tag to be considered
     * @param selected_nodes mask to the nodes to be considered (optional)
     * @return map of tags to nodes
     */
    std::map<bsg10xTag, std::vector<sgNodeID_t>> get_tag_nodes(uint32_t min_nodes = 2,
                                                               const std::vector<bool> &selected_nodes = {});

    /** @brief Returns a collection of pairs of nodeIDs that share more than min_shared tags between them
     * Returns a list of pairs of nodes that are neighbours according to 10x reads.
     * This is th eone that makes the reflective thing
     * @param min_shared
     * @param selected_nodes
     * @return
     */
    std::vector<std::pair<sgNodeID_t , sgNodeID_t >> get_tag_neighbour_nodes(uint32_t min_shared,const std::vector<bool> & selected_nodes={});

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
     */
    void compute_all_tag_neighbours(int min_size,float min_score);
    void compute_all_tag_neighbours2(int min_size,float min_score, int min_mapped_reads_per_tag=2);

    void write_tag_neighbours(std::string filename);
    void read_tag_neighbours(std::string filename);

    const SequenceDistanceGraph & sg;
    const WorkSpace &ws;

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
     *  Collection of neightbouring nodes
     *  Completed using compute_all_tag_neighbours()
     *  tag_neighbours[nodeID] = [collection of 10 determined neighbours (see TagNeighbour) ...]
     */
    std::vector<std::vector<TagNeighbour>> tag_neighbours; //not persisted yet!

    static const bsgVersion_t min_compat;
};
