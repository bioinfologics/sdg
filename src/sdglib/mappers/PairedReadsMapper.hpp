//
// Created by Bernardo Clavijo (EI) on 12/05/2018.
//

#ifndef BSG_PAIREDREADSMAPPER_HPP
#define BSG_PAIREDREADSMAPPER_HPP

#include <map>
#include <fstream>

#include "sdglib/types/MappingTypes.hpp"
#include "sdglib/factories/KMerIDXFactory.hpp"
#include "sdglib/readers/SequenceGraphReader.hpp"
#include <sdglib/datastores/PairedReadsDatastore.hpp>
#include <sdglib/indexers/UniqueKmerIndex.hpp>

class UniqueKmerIndex;
class Unique63merIndex;
class PairedReadConnectivityDetail; //Forward declaration
class WorkSpace;

/**
 * @brief A mapper for linked reads from a PairedReadsDatastore.
 *
 * Supports partial remapping of unmapped reads or of a selection list.
 */

class PairedReadsMapper {

public:
    PairedReadsMapper(WorkSpace &_ws, PairedReadsDatastore &_datastore);
    void write(std::ofstream & output_file);
    void read(std::ifstream & input_file);
    /** @brief Maps each read in the data-store to the nodes using unique kmers form the graph
     *
     * Maps each read in the data-store to the nodes using unique kmers form the graph.
     * Reads are mapped only if the placing in the graph is unique, if one or more kmers of a read
     * maps to a different node the mapping for that read is discarded.
     * Results are stored in this.read_to_node and this.reads_in_node
     *
     * If a set of reads ids is passed the mapper will only remap those reads.
     *
     * @param reads_to_remap
     */
    void map_reads(std::unordered_set<uint64_t> const &  reads_to_remap={});
    /**
     * Clears this.read_to_node and this.reads_in_node collections and reruns map_reads()
     */
    void remap_all_reads();

    /** @brief same as map_reads() but using k=63
     *
     * @param reads_to_remap
     */
    void map_reads63(std::unordered_set<uint64_t> const &  reads_to_remap={});

    /** @brief Same as remap_all_reads but using map_reads63
     *
     */
    void remap_all_reads63();

    /** @brief Discards mappings of nodes marked as deleted.
     *
     * Discards mappings of nodes marked as deleted by changing to 0 the mappings in read_to_node[read_id] and clearing the reads_in_node[deletedNode] vector
     */
    void remove_obsolete_mappings();

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
    void print_status();

    PairedReadsMapper operator=(const PairedReadsMapper &other);

    /** @brief Returns a collection of read ids that have both ends mapped to the nodeID.
     *
     * @param nodeID id of the node to get the read ids from
     * @return vector of read ids where both ends are mapped to the nodeID
     */
    std::vector<uint64_t> get_node_readpairs_ids(sgNodeID_t nodeID);

    const SequenceDistanceGraph & sg;
    const std::shared_ptr<UniqueKmerIndex> kmer_to_graphposition;
    const std::shared_ptr<Unique63merIndex> k63mer_to_graphposition;
    const PairedReadsDatastore & datastore;
    /**
     * reads_in_node[i] contains all the read ids from the data-store mapped to node i
     */
    std::vector<std::vector<ReadMapping>> reads_in_node;

    /**
     * read_to_node[i] containes the node where read i was mapped
     */
    std::vector<sgNodeID_t> read_to_node; //id of the main node if mapped, set to 0 to remap on next process
    //TODO: reading and writing this would simplify things??
    /**
     * read_direction_in_node[i] has the direction with wich read i was mapped in the corresponding mapping (in read_to_node[i])
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

    static const bsgVersion_t min_compat;

};

/**
 * @brief Analysis of all reads connecting two particular nodes.
 */

class PairedReadConnectivityDetail {
public:
    PairedReadConnectivityDetail(){};
    PairedReadConnectivityDetail(const PairedReadsMapper & prm, sgNodeID_t source, sgNodeID_t dest);
    PairedReadConnectivityDetail& operator+=(const PairedReadConnectivityDetail& rhs){
        this->orientation_paircount[0] += rhs.orientation_paircount[0];
        this->orientation_paircount[1] += rhs.orientation_paircount[1];
        this->orientation_paircount[2] += rhs.orientation_paircount[2];
        this->orientation_paircount[3] += rhs.orientation_paircount[3];
        return *this;
    }

    uint64_t orientation_paircount[4]={0,0,0,0};
};

#endif //BSG_PAIREDREADSMAPPER_HPP
