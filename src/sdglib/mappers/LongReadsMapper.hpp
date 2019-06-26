//
// Created by Luis Yanes (EI) on 12/02/2018.
//

#ifndef BSG_LONGREADMAPPER_H
#define BSG_LONGREADMAPPER_H


#include <iostream>
#include <memory>
#include <vector>
#include <sdglib/types/GenericTypes.hpp>
#include <sdglib/indexers/NKmerIndex.hpp>
#include <sdglib/indexers/SatKmerIndex.hpp>
#include <sdglib/types/MappingTypes.hpp>
#include <sdglib/utilities/hashing_helpers.hpp>
#include <sdglib/Version.hpp>
#include <sdglib/datastores/ReadSequenceBuffer.hpp>


class WorkSpace;
class LinkedReadsMapper;
struct ReadPathParams {
    int default_overlap_distance = 199;
    float path_distance_multiplier = 1.5;
    int min_path_distance = 400;
    int min_negative_path_distance = 500;
    int max_path_nodes = 40;
    int max_distinct_paths = 1000;
    int path_mapping_k = 15;
    int max_mapping_path_size = 1000;
    int min_mapping_chain = 2;
    int max_mapping_jump = 100;
    int max_mapping_delta_change = 30;
    float kmer_filter_multiplier = 2;
};

struct ReadCacheItem {
    uint64_t id;
    std::string seq;
    ReadCacheItem(uint64_t i, std::string &sequence): id(i), seq(sequence) {}

    friend std::ostream& operator<<(std::ostream& os, const ReadCacheItem& read) {
        os << read.id;
        return os;
    }

    bool operator<(const ReadCacheItem &o) const {
        return id<o.id;
    }

};

struct match_band{
    bool dir;
    int32_t min_offset;
    int32_t max_offset;
    int32_t len;
    int32_t score;

    match_band(bool d, int32_t s, int32_t e, int32_t l, int32_t sc) :
    dir(d),
    min_offset(s),
    max_offset(e),
    len(l),
    score(sc){}
};

class HaplotypeScore {
public:
    float score;
    std::vector<sgNodeID_t> haplotype_nodes;

    bool operator<(const HaplotypeScore &other) const {
        return score<other.score;
    }
};

class LongReadsMapper;
class LongReadsDatastore;

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
    std::vector<uint64_t> rkmers;
    std::vector<uint64_t> nkmers;
    std::vector<uint8_t> coverage;
    std::vector<std::vector<bool>> kmer_hits_by_node;
public:
    LongReadHaplotypeMappingsFilter (const LongReadsMapper & _lorm, const LinkedReadsMapper & _lirm);
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
    const LongReadsMapper & lorm;
    const LinkedReadsMapper & lirm;
    ReadSequenceBuffer * lrbsgp;
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
class LongReadsMapper {
    NKmerIndex assembly_kmers;
    SatKmerIndex sat_assembly_kmers;
public:

    const SequenceDistanceGraph & sg;

    // TODO: Describe exactly what each of these do!!
    uint8_t k=15;
    int min_size=1000;
    int min_chain=50;
    int max_jump=500;
    int max_delta_change=60;

    LongReadsMapper(const WorkSpace &_ws, const LongReadsDatastore &ds, uint8_t k=15, bool sat_index=false);
    LongReadsMapper(const SequenceDistanceGraph &_sdg, const LongReadsDatastore &ds, uint8_t k=15, bool sat_index=false);

    LongReadsMapper& operator=(const LongReadsMapper &other);
    LongReadsMapper(const LongReadsDatastore &ds, const LongReadsMapper &o);

    void print_status() const;

    /** @brief Getter for the defined datastore
     *
     * @return
     */
    const LongReadsDatastore& getLongReadsDatastore() {return datastore;}

    /** @brief Getter for the defined SequenceGraph
     *
     * @return
     */
    const SequenceDistanceGraph& getSequenceGraph() {return sg;}

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
     * @param filter_limit
     * @param readIDs
     * @param detailed_log
     */
    void map_reads(int filter_limit = 200, const std::unordered_set<uint32_t> &readIDs = {},
                   std::string detailed_log = "");

    /**
     * This function maps any sequence to the graph, index needs to be already updated!
     * WARNING: this is slow, not meant for high throughput
     * @param query_sequence_ptr
     * @return
     */
    std::vector<LongReadMapping> map_sequence(const char * query_sequence_ptr, sgNodeID_t seq_id=0);
    //void map_reads(std::string detailed_log){map_reads({},detailed_log);};

    void read(std::string filename);

    void read(std::ifstream &inf);

    void write(std::string filename);

    void write(std::ofstream &output_file);

    void write_filtered_mappings(std::string filename);

    void read_filtered_mappings(std::string filename);

    void write_read_paths(std::string filename);

    void read_read_paths(std::string filename);

    /**
     * Updates reads_in_nodes
     */
    void update_indexes();

    /**
     * Updates the assembly_kmers index with the kmers of the current graph with frequency less than 200
     */
    void update_graph_index(int filter_limit=200, bool verbose=true);

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
    void filter_mappings_with_linked_reads(const LinkedReadsMapper &lrm, uint32_t min_size=10000, float min_tnscore=0.03, uint64_t first_id=0, uint64_t last_id=0);

    void filter_mappings_by_size_and_id(int64_t size, float id);
    /**
     * This goes read by read, and transforms the filtered mappings into path-mappings, by:
     * 1) chaining coherent successive mappings
     * 2) removing secondary mappings
     * 3) Pathing between mappings as far as the read supports a single path
     *
     */
//    void path_filtered_mappings();

    //void dump_path_mappings();

    //void load_path_mappings();

    std::vector<LongReadMapping>
    filter_and_chain_matches_by_offset_group(std::vector<LongReadMapping> &matches, bool verbose=false);

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
    std::vector<LongReadMapping> improve_read_filtered_mappings(uint32_t rid, bool correct_on_ws=false);

    void improve_filtered_mappings() {
#pragma omp parallel for
        for (uint32_t rid = 0; rid < filtered_read_mappings.size(); ++rid) {
            improve_read_filtered_mappings(rid,true);
        }
        update_indexes();
    }

    std::vector<ReadCacheItem>create_read_paths(const std::vector<sgNodeID_t> &backbone, const ReadPathParams &read_path_params);

    std::vector<sgNodeID_t> create_read_path(uint32_t rid, const ReadPathParams &read_path_params, bool verbose=false, const std::string& read_seq="");

//    std::vector<sgNodeID_t> create_read_path_fast(uint32_t rid, bool verbose=false, const std::string read_seq="");
    /**
     * This updates the filtered mappings by taking the elements from the path mapping that have the right size and neighbourhood conditions
     */
    //void update_filtered_mappings_from_paths(const LinkedReadsMapper &lrm,uint32_t min_size=10000, float min_tnscore=0.03);

    //void update_filtered_mappings_from_paths(uint32_t min_size=10000, float min_tnscore=0.03);

    const LongReadsDatastore &datastore;
    /**
     * Flat list of mappings from the reads
     */
    std::vector<LongReadMapping> mappings;
    std::vector<int64_t> first_mapping; //index to the first mapping of the read. If no mappings, -1.

    /**
     * This structure holds "all paths" between consecutive backbone anchors
     *
     * (Maybe in canonical from->to orientation?)
     */
    std::unordered_map<std::pair<sgNodeID_t,sgNodeID_t>, std::vector<SequenceGraphPath>> all_paths_between;

    std::vector < std::vector<LongReadMapping> > filtered_read_mappings;

    std::vector<std::vector<sgNodeID_t>> read_paths;

    /**
     * Stores an index of all reads that map to a node.
     * This index can be queried restrict search of particular mappings.
     *
     */
    std::vector<std::vector<uint64_t>> reads_in_node;


    static const sdgVersion_t min_compat;

    bool sat_kmer_index = false;
};


#endif //BSG_LONGREADMAPPER_H
