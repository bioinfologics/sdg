//
// Created by Bernardo Clavijo (EI) on 2019-07-17.
//

#pragma once

#include <sdglib/graph/DistanceGraph.hpp>
#include <sdglib/workspace/WorkSpace.hpp>

class HaplotypeScore {
public:
    float score;
    std::vector<sgNodeID_t> haplotype_nodes;

    bool operator<(const HaplotypeScore &other) const {
        return score<other.score;
    }
};

class LongReadsMapper;

/**
 * This class groups all methods to filter long read mappings to haplotype solutions within a long read
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
    /**
     * From an initiated read instance creates all unique sets of haplotipic neighbours using the linked reads neighbour matrix
     * The result is stores in haplotype_scores in the class
     * @param min_tn min score for 2 nodes to be considered neighbours
     */
    void generate_haplotypes_from_linkedreads(float min_tn=0.05);
    /**
     * Scores each of the haplotype scores according to how much of the read the haplotype covers, the results are added to the score attibute inside each element of haplotype_scores
     * @param weight
     */
    void score_coverage(float weight);
    /**
     *
     * @param weight
     * @param k
     * @param win_size
     * @param win_step
     */
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

/**
 * @brief Creates a new DistanceGraph with linkage information from the workspace over the nodes of an input DistanceGraph
 *
 */
class LinkageMaker {
public:
    explicit LinkageMaker(const DistanceGraph &_dg):dg(_dg){deselect_all();};
    explicit LinkageMaker(const SequenceDistanceGraph &_dg):LinkageMaker(static_cast<const DistanceGraph &>(_dg)){};

    void check_selected_nodes();
    void deselect_all();
    void select_all();
    void report_selection();
    void select_by_size(uint64_t min_size,uint64_t max_size=0);

    /**
     * Creates a distance graph following the topology of the starting graph
     *
     * A->b->c->D->e->f->g->H
     *
     * Uppercase letters are selecetd nodes
     *
     * if called wwith radius=10 A-D y D-H are connected
     *
     * if called with radius=3 A is connected to D, but D is not connected to H
     *
     * if called with radius=2, nothing gets connected
     *
     *
     * @param radius reach of the connections
     * @return
     */
    DistanceGraph make_topology_linkage(int radius);

    /**
     * Create a distance graph using the available paired read mappers
     * @param min_reads minimum amount of reads to create a link
     * @return
     */
    DistanceGraph make_paired_linkage(int min_reads);

    DistanceGraph make_paired_linkage_pe(int min_reads);

    //Linkage creation methods (work on selected nodes)
    std::map<std::pair<sgNodeID_t, sgNodeID_t>, uint64_t> shared_read_paths(int min_shared, std::vector<size_t> libraries, bool r1rev, bool r2rev);

    //supporting methods
    std::vector<Link> mappings_to_multilinkage(const std::vector<LongReadMapping> &lorm_mappings, uint32_t read_size, int32_t unmapped_end=1000);

    DistanceGraph make_longreads_multilinkage(const LongReadsMapper &lorm, uint64_t min_map_size=1000, float min_map_id=.1, bool real_read_size=true, int32_t unmapped_end=1000);

//    DistanceGraph make_longreads_multilinkage(const std::string &datastore_name, uint64_t min_map_size=1000, float min_map_id=.1, bool real_read_size=true, int32_t unmapped_end=1000);
    /**
     * Creates a multilinkage digraph using the Nano10X pipeline
     * @param lorm long read mapper
     * @param lrm link read mapper (needs to have the neighbours calculated, see linkedReadsMapper.compute_all_tag_neighbours)
     * @param min_size min node size to be considered
     * @param min_tnscore min neighbour tag score to be considered
     * @param real_read_size
     * @param unmapped_end
     * @return Multi linkage DistanceGraph
     */
    DistanceGraph make_long10x_multilinkage(const LongReadsMapper &lorm, const LinkedReadsMapper &lrm, uint32_t min_size,  float min_tnscore, bool real_read_size=true, int32_t unmapped_end=1000);

    static std::vector<std::vector<LongReadMapping>> filter_mappings_with_linked_reads(const LongReadsMapper &lorm, const LinkedReadsMapper &lrm, uint32_t min_size,  float min_tnscore);

    DistanceGraph make_long10x_multilinkage(const std::string &lorm_name, const std::string &lrm_name, uint32_t min_size,  float min_tnscore, bool real_read_size=true, int32_t unmapped_end=1000);

    DistanceGraph make_paired10x_multilinkage(const PairedReadsMapper &prm, const LinkedReadsMapper &lirm, float min_tnscore=0.2, bool fr=false, uint64_t read_offset=0, int64_t min_dist_size=-1);


    const DistanceGraph &dg;
    std::vector<bool> selected_nodes;
};

