//
// Created by Bernardo Clavijo (EI) on 2019-07-12.
//

#pragma once

#include <sdglib/indexers/NKmerIndex.hpp>
#include <sdglib/indexers/SatKmerIndex.hpp>
#include <sdglib/types/MappingTypes.hpp>

/**
 * This class implements a simple sequence mapper to query a DG with single sequences.
 * On instanciation, this class builds an index.
 * It provides a map function that returns a SequenceMappingsView for sequence.
 */
class SequenceMapper {
public:
    const DistanceGraph & sg;

    // TODO: Describe exactly what each of these do!!
    uint8_t k=15;
    int min_size=1000;
    int min_chain=50;
    int max_jump=500;
    int max_delta_change=60;

    explicit SequenceMapper(const DistanceGraph &_dg,uint8_t _k=15, int _min_size=1000, int _min_chain=50, int _max_jump=500, int _max_delta_change = 60);
    explicit SequenceMapper(const SequenceDistanceGraph &_dg,uint8_t _k=15, int _min_size=1000, int _min_chain=50, int _max_jump=500, int _max_delta_change = 60):
        SequenceMapper(static_cast<const DistanceGraph &>(_dg),_k,_min_size,_min_chain,_max_jump,_max_delta_change){};//This is just for SWIG's benefit

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
     * This function maps any sequence to the graph, index needs to be already updated!
     * WARNING: this is slow, not meant for high throughput
     * @param query_sequence_ptr
     * @return
     */
    std::vector<LongReadMapping> map_sequence(const char * query_sequence_ptr, sgNodeID_t seq_id=0);

    /**
     * Updates the assembly_kmers index with the kmers of the current graph with frequency less than 200
     */
    void update_graph_index(int filter_limit=200, bool verbose=true);

    /**
     * This goes read by read, and filters the mappings by finding a set of linked nodes that maximises 1-cov of the read
     *
     * Unfiltered mappings read from mappings and results stored in filtered_read_mappings, which is cleared.
     *
     * @param lrm a LinkedReadMapper with mapped reads, over the same graph this mapper has mapped Long Reads.
     * @param min_size minimum size of the read to filter mappings.
     * @param min_tnscore minimum neighbour score on linked reads
     * @param first_id
     * @param last_id
     */
    //void filter_mappings_by_size_and_id(int64_t size, float id);

    /**
     * This goes read by read, and transforms the filtered mappings into path-mappings, by:
     * 1) chaining coherent successive mappings
     * 2) removing secondary mappings
     * 3) Pathing between mappings as far as the read supports a single path
     *
     */
    //std::vector<LongReadMapping>
    //filter_and_chain_matches_by_offset_group(std::vector<LongReadMapping> &matches, bool verbose=false);

    /**
     * Eliminates matches that are contained within another bigger better match (more span and more score)
     * @param matches set of matched to de-shadow, usually all the matches within a read
     * @param verbose
     * @return vector of de-shadowed LongReadMappings
     */
    std::vector<LongReadMapping>
    remove_shadowed_matches(std::vector<LongReadMapping> &matches, bool verbose=false);

    /**
     * Performs the de-shadowing process for each read
     * @param rid
     * @param correct_on_ws
     * @return
     */
    //std::vector<LongReadMapping> improve_read_filtered_mappings(uint32_t rid, bool correct_on_ws=false);


private:
    /**
     * Populates the matches container with the matches between all kmers from one read and the *saturated* index of the graph.
     * The index is a filtered set of kmers from the graph constructed using update_graph_index() or similar.
     * A match is a perfect Kmer match between the read and the graph index
     *
     * matches are stored as kmer_index_in_read->kmer_match->match_description:
     *  matches[kmer_index_in_read][kmer_match].first = node_id for a match of that kmer in the graph, sign indicates orientation
     *  matches[kmer_index_in_read][kmer_match].second = offset for a match of that kmer for corresponding node_id (.first), the offsets are always calculated from the begining of the read
     *
     *  TODO: Change match type from in32_t to sgNodeId in this function !!!
     *
     * @param matches structure to store kmer mappings
     * @param read_kmers contains read kmers with orientations
     */
    void get_sat_kmer_matches(std::vector<std::vector<std::pair<int32_t, int32_t>>> &matches, std::vector<std::pair<bool, uint64_t>> &read_kmers);

    /**
     * Populates the matches container with the matches between all kmers from one read and the index of the graph.
     * The index is a filtered set of kmers from the graph constructed using update_graph_index() or similar.
     * A match is a perfect Kmer match between the read and the graph index
     *
     * matches are stored as kmer_index_in_read->kmer_match->match_description:
     *  matches[kmer_index_in_read][kmer_match].first = node_id for a match of that kmer in the graph, sign indicates orientation
     *  matches[kmer_index_in_read][kmer_match].second = offset for a match of that kmer for corresponding node_id (.first), the offsets are always calculated from the begining of the read
     *
     *  TODO: Change match type from in32_t to sgNodeId in this function !!!
     *
     * @param matches structure to store kmer mappings
     * @param read_kmers contains read kmers with orientations
     */
    void get_all_kmer_matches(std::vector<std::vector<std::pair<int32_t, int32_t>>> & matches, std::vector<std::pair<bool, uint64_t>> & read_kmers);

    /**
     * Using the populated matches collection (see het_all_kmer_matches()) returns a collection of nodeIDs that have more than 50 matches with that particular read.
     * TODO: parameter to ckeck, the 50 in the filter, why is 50?.
     * @param matches matches vector
     * @param read_kmers_size number of kmers in the read
     * @return set of nodeIDs of nodes with more than 50 matches to the read
     */
    void
    count_candidates(std::vector<unsigned char> &candidate_counts,
                     std::vector<std::vector<std::pair<int32_t, int32_t>>> &matches,
                     uint32_t read_kmers_size);

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
    std::vector<LongReadMapping> alignment_blocks(uint32_t readID,
                                                  std::vector<std::vector<std::pair<int32_t, int32_t>>> &matches,
                                                  uint32_t read_kmers_size, const std::vector<unsigned char> &candidate_counts);

    /**
     * Given a list of blocks the filter will discard overlapping blocks keeping those with the max span and score combination
     *
     * @param blocks Maping blocks produced by alignment_blocks()
     * @param matches Matches produced by get_all_kmer_matches()
     * @param read_kmers_size Number of kmers in the read
     * @return Filtered, non overlapping collection of mapping blocks (LongReadMapping)
     */
    std::vector<LongReadMapping> filter_blocks(std::vector<LongReadMapping> & blocks, std::vector<std::vector<std::pair<int32_t, int32_t>>> & matches,  uint32_t read_kmers_size);


    NKmerIndex assembly_kmers;
    SatKmerIndex sat_assembly_kmers;

    std::vector<LongReadMapping> filtered_mappings;
    bool sat_kmer_index = false;
};
