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
#include <memory>

class HaplotypeScore {
public:
    float score;
    std::vector<sgNodeID_t> haplotype_nodes;

    bool operator<(const HaplotypeScore &other) const {
        return score<other.score;
    }
};

class LongReadMapper;

/**
 * This class groups all methods to filter long read mappings to  haplotype solutions within a long read
 * The
 * 1) set_read -> inits the problem, gets read sequence, mappings. cleans up alignments if needed.
 * 2) generate_haplotypes_from_linked_reads -> creates the possible haplotypes, all with score 0
 * 3) score_* -> these funcions add up scores in range(0,1), modulated by their weight parameter
 * 4) sort
 * 5) filter_winner_haplotype -> populates filtered_alingments with the winning haplotype alignments
 */
class LongReadHaplotypeMappingsFilter {
public:
    LongReadHaplotypeMappingsFilter (const LongReadMapper & _lorm, const LinkedReadMapper & _lirm);
    ~LongReadHaplotypeMappingsFilter(){
        delete(lrbsgp);
    }
    void set_read(uint64_t read_id);
    void generate_haplotypes_from_linkedreads(float min_tn=0.05);
    void score_coverage(float weight);
    void score_window_winners(float weight, int k=15, int win_size=500, int win_step=250);

    /**
     * This method runs all the scoring, then puts all haplotypes that pass criteria and their scores in a sorted vector
     *
     * @param read_id
     * @param coverage_weight
     * @param winners_weight
     * @param min_tn
     * @return
     */
    void rank(uint64_t read_id, float coverage_weight, float winners_weight, float min_tn=0.05);

    uint64_t read_id;
    std::string read_seq;
    std::vector<LongReadMapping> mappings;
    std::vector<HaplotypeScore> haplotype_scores;
    const LongReadMapper & lorm;
    const LinkedReadMapper & lirm;
    BufferedSequenceGetter * lrbsgp;
    std::vector<sgNodeID_t> nodeset;


};

enum MappingFilterResult {Success, TooShort, NoMappings, NoReadSets, LowCoverage};

/**
 * Long read mapping to SequenceGraph, computation and storage of the raw alignments and filtered alingments.
 *
 * this->mappings is filled via small k-mers to multi-position index and a chain search.
 * this->filtered_read_mappings is filled by calling one of the filter_mappings_* methods, which can use extra data.
 * the reads_in_node index is populated by update_indexes() from this->filtered_read_mappings data.
 */
class LongReadMapper {
    NKmerIndex assembly_kmers;

public:

    const SequenceGraph & sg;

    uint8_t k=15;
    int min_size=1000;
    int min_chain=50;
    int max_jump=500;
    int max_delta_change=60;

    LongReadMapper(SequenceGraph &sg, LongReadsDatastore &ds, uint8_t k=15);
    ~LongReadMapper();

    LongReadMapper operator=(const LongReadMapper &other);

    /** @brief Getter for the defined datastore
     *
     * @return
     */
    LongReadsDatastore& getLongReadsDatastore() {return datastore;}

    /** @brief Getter for the defined SequenceGraph
     *
     * @return
     */
    const SequenceGraph& getSequenceGraph() {return sg;}

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
     * Copies all mappings to read from raw mappings into a vector and returns it.
     * @param read_id
     * @return
     */
    std::vector<LongReadMapping> get_raw_mappings_from_read(uint64_t read_id) const;

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
     * Updates reads_in_nodes
     */
    void update_indexes();

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
     * @param first_id
     * @param last_id
     */
    void filter_mappings_with_linked_reads(const LinkedReadMapper &lrm, uint32_t min_size=10000, float min_tnscore=0.03, uint64_t first_id=0, uint64_t last_id=0);

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
     * Flat list of mappings from the reads
     */
    std::vector<LongReadMapping> mappings;
    std::vector<int64_t> first_mapping; //index to the first mapping of the read. If no mappings, -1.

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
