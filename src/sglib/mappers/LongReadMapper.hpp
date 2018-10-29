//
// Created by Luis Yanes (EI) on 12/02/2018.
//

#ifndef BSG_LONGREADMAPPER_H
#define BSG_LONGREADMAPPER_H


#include <iostream>
#include <sglib/factories/KMerIDXFactory.hpp>
#include <sglib/datastores/LongReadsDatastore.hpp>
#include <sglib/graph/SequenceGraph.hpp>
#include <sglib/types/MappingTypes.hpp>
#include <sglib/indexers/NKmerIndex.hpp>
#include "LinkedReadMapper.hpp"


enum MappingFilterResult {Success, TooShort, NoMappings, NoReadSets, LowCoverage};
/**
 * Long read mapping to the graph, this class manages storage and computation of the alignments.
 */
class LongReadMapper {

    const SequenceGraph & sg;

    uint8_t k=15;
    int min_size=1000;
    int min_chain=50;
    int max_jump=500;
    int max_delta_change=60;

    NKmerIndex assembly_kmers;

    void update_indexes();

public:

    LongReadMapper(SequenceGraph &sg, LongReadsDatastore &ds, uint8_t k=15);
    ~LongReadMapper();

    LongReadMapper operator=(const LongReadMapper &other);

    /** @brief Getter for the defined datastore
     *
     * @return
     */
    LongReadsDatastore& getLongReadsDatastore() {return datastore;}

    /**
     * Sets mapping parameters
     * TODO: explain each parameter
     * @param _k kmer size for mapping
     * @param _min_size
     * @param _min_chain
     * @param _max_jump
     * @param _max_delta_change
     */
    void set_params(uint8_t _k=15, int _min_size=1000, int _min_chain=50, int _max_jump=500, int _max_delta_change = 60){
        k=_k;
        min_size=_min_size;
        min_chain=_min_chain;
        max_jump=_max_jump;
        max_delta_change=_max_delta_change;
    }

    /**
     * Populates the matches container with the matches between all kmers from one read and the index of the graph.
     * The index is a filtered set of kmers from the graph constructed using update_graph_index() or similar.
     * A match is a perfect Kmer match between the read and the graph index
     *
     * matches are stored as read_kmer->kmer_match->match_description:
     *  matches[read_kmer][kmer_match].first = node_id for a match of that kmer in the graph, sign indicates orientation
     *  matches[read_kmer][kmer_match].second = offset for a match of that kmer for corresponding node_id (.first), the offsets are always calculated from the begining of the read
     *
     *  TODO: Change matche type from in32_t to sgNodeId in this function !!!
     *
     * @param vetor to store the mappings
     * @param vector containing read kmers with the orientations
     */
    void get_all_kmer_matches(std::vector<std::vector<std::pair<int32_t, int32_t>>> & matches, std::vector<std::pair<bool, uint64_t>> & read_kmers);

    /**
     * Using the populated matches collection (see het_all_kmer_matches()) returns a collection of nodeIDs that have more than 50 matches with that particular read.
     * TODO: parameter to ckeck, the 50 in the filter, why is 50?.
     * @param matches matches vector
     * @param read_kmers_size number of kmers in the read
     * @return set of nodeIDs of nodes with more than 50 matches to the read
     */
    std::set<sgNodeID_t> window_candidates(std::vector<std::vector<std::pair<int32_t, int32_t>>> & matches, uint32_t read_kmers_size);

    /**
     * This function will check for continuous blocks of alignments between the read and the window candidates produced by window_candidates() using the matches prudiced by get_allKmer_matches()
     * The result is a vector af LongReadMappings where each element describes a mapping of a block of the read to a block of a node (see struct LongReadMapping)
     *
     * Each LongReadMapping describes a chain of matching kmers in the read and the node where the matches advance at the same pace (within a jumping distance allowance)
     *
     * @param readID id of the analyzed read
     * @param matches matches collection produced by get_allKmer_matches()
     * @param read_kmers_size Number of kmers in the read
     * @param candidates collection of candidate matches produced by window_candidates()
     * @return collection of readtonode mappings (LongReadMapping)
     */
    std::vector<LongReadMapping> alignment_blocks(uint32_t readID, std::vector<std::vector<std::pair<int32_t, int32_t>>> & matches,  uint32_t read_kmers_size, std::set<sgNodeID_t> &candidates);

    /**
     * Given a list of blocks the filter will discard overlapping blocks keeping those with the max span and score combination
     *
     * @param blocks Maping blocks produced by alignment_blocks()
     * @param matches Matches produced by get_all_kmer_matches()
     * @param read_kmers_size Number of kmers in the read
     * @return Filtered, non overlapping collection of mapping blocks (LongReadMapping)
     */
    std::vector<LongReadMapping> filter_blocks(std::vector<LongReadMapping> & blocks, std::vector<std::vector<std::pair<int32_t, int32_t>>> & matches,  uint32_t read_kmers_size);

    /**
     * Function to map a read to the graph in 4 steps using the methods in this class
     *
     * //========== 1. Get read sequence, kmerise, get all matches ==========
     * //========== 2. Find match candidates in fixed windows ==========
     * //========== 3. Create alignment blocks from candidates ==========
     * //========== 4. Construct mapping path ==========
     *
     * Results are stored in the mappings collection of this object
     *
     * @param readIDs
     * @param detailed_log
     */
    void map_reads(std::unordered_set<uint32_t> readIDs = {},std::string detailed_log="");

    void map_reads(std::string detailed_log){map_reads({},detailed_log);};

    void read(std::string filename);

    void read(std::ifstream &inf);

    void write(std::string filename);

    void write(std::ofstream &output_file);

    void write_filtered_mappings(std::string filename);

    void read_filtered_mappings(std::string filename);

    /**
     * Updates the assembly_kmers index with the kmers of the current graph with frequency less than 200
     */
    void update_graph_index();

    /**
     * This goes read by read, and filters the mappings by finding a set of linked nodes that maximises 1-cov of the read
     *
     * Unfiltered mappings read from this->mappings and results stored in this->filtered_read_mappings, which is cleared.
     *
     * @param lrm a LinkedReadMapper with mapped reads, over the same graph this mapper has mapped Long Reads.
     * @param min_size minimum size of the read to filter mappings.
     * @param min_tnscore minimum neighbour score on linked reads
     */
    void filter_mappings_with_linked_reads(const LinkedReadMapper &lrm, uint32_t min_size=10000, float min_tnscore=0.03);

    /**
     * Single-read nano10x filtering. Return status
     * @param lrm
     * @param lrbsg
     * @param min_size
     * @param min_tnscore
     * @param readID
     * @param offset_hint
     * @return
     */
    MappingFilterResult filter_mappings_with_linked_reads(const LinkedReadMapper &lrm, BufferedSequenceGetter &lrbsg, uint32_t min_size, float min_tnscore, uint64_t readID, uint64_t offset_hint=0);

    LongReadsDatastore datastore;
    /**
     * This public member stores a flat list of mappings from the reads, it is accessed using the mappings_in_node index
     * or the read_to_mappings index.
     */
    std::vector<LongReadMapping> mappings;
    std::vector < std::vector<LongReadMapping> > filtered_read_mappings;

    /**
     * Stores an index of all reads that map to a node.
     * This index can be queried restrict search of particular mappings.
     *
     */
    std::vector<std::vector<uint64_t>> reads_in_node;




    static const bsgVersion_t min_compat;

};


#endif //BSG_LONGREADMAPPER_H
