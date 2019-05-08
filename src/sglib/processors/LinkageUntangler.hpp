//
// Created by Bernardo Clavijo (EI) on 28/05/2018.
//

#ifndef BSG_LINKAGEUNTANGLER_HPP
#define BSG_LINKAGEUNTANGLER_HPP


#include <sglib/workspace/WorkSpace.hpp>
#include <sglib/graph/LinkageDiGraph.hpp>

/**
 * @brief Generates and manipulates links from 10x, paired and long reads
 *
 * Uses a workspace containing mapped reads and KCI indexes.
 *
 * This implements a KISS base for graph simplification/scaffolding in three steps:
 *
 * 1) choose appropriate nodes
 * 2) link them
 * 3) improve the linkage by locally analysing the choices from (1) and (2)
 *
 * It also provides a number of functions to modify the graph according to linkage, which should not be here.
 *
 */
class LinkageUntangler {
public:

    explicit LinkageUntangler(WorkSpace & _ws): ws(_ws) { clear_node_selection();};

    //==== Methods currently in use =====

    //Node selection methods
    void clear_node_selection();
    void report_node_selection();
    void select_nodes_by_size_and_ci( uint64_t min_size, float min_ci, float max_ci);
    std::set<std::pair<sgNodeID_t, sgNodeID_t >> get_HSPNPs(uint64_t min_size, float min_ci, float max_ci);

    /**
    * Anchors are chosen by having all their strongly connected neighbours (i.e. those with min_links connections)
    * coherently connected among them in order with at least min_transitive_links.
    * @param multi_ldg
    * @param min_links
    * @param min_transitive_links
    */
    void select_multi_linkage_linear_anchors(const LinkageDiGraph & multi_ldg, int min_links=3, int min_transitive_links=2);

    //Multi-Linkage creation methods: multiple evidence-supported links between any nodes

    LinkageDiGraph make_topology_linkage(int radius);
    LinkageDiGraph make_longRead_multilinkage(const LongReadMapper &lorm,bool real_read_size=true, int32_t unmapped_end=1000);
    LinkageDiGraph make_paired10x_multilinkage(const PairedReadMapper &prm, const LinkedReadMapper &lirm, float min_tnscore=0.2, bool fr=false, uint64_t read_offset=0);


    //supporting methods
    std::vector<Link> mappings_to_multilinkage(const std::vector<LongReadMapping> &lorm_mappings, uint32_t read_size, int32_t unmapped_end=1000);



    //Linkage creation methods: aggregated links between selected nodes (anchors)

    LinkageDiGraph make_nextselected_linkage(const LinkageDiGraph & multi_ldg, int min_links=3);

    LinkageDiGraph make_and_simplify_anchor_linkage(const LinkageDiGraph &multi_ldg, int anchor_link=5, int transitive_links=3);



    //==== DEPRECATED methods, do NOT use in new projects =====


    void select_nodes_by_HSPNPs(uint64_t min_size, float min_ci, float max_ci);



    //Linkage creation methods (work on selected nodes)
    std::map<std::pair<sgNodeID_t, sgNodeID_t>, uint64_t> shared_read_paths(int min_shared, std::vector<size_t> libraries, bool r1rev, bool r2rev);

    LinkageDiGraph make_paired_linkage(int min_reads);
    LinkageDiGraph make_paired_linkage_pe(int min_reads);
    LinkageDiGraph make_tag_linkage(int min_tags, bool use_kmer_paths=false);

    LinkageDiGraph make_and_simplify_linkage(int min_shared_tags);

    //Linkage improving/filtering methods
    LinkageDiGraph filter_linkage_to_hspnp_duos( uint64_t min_size, float min_ci, float max_ci, const LinkageDiGraph & ldg);


    //Graph untangling/modification/local assembly methods
    void expand_trivial_repeats(const LinkageDiGraph &);
    void expand_linear_regions(const LinkageDiGraph &);
    void expand_linear_regions_skating(const LinkageDiGraph &, int max_lines=0);
    void linear_regions_tag_local_assembly(const LinkageDiGraph & ldg, uint8_t k, int min_cvg, int max_lines, uint64_t min_nodes, uint64_t min_total_size, bool count_tag_cvg=false);

    //Problem localisation methods
    void fill_linkage_line(std::vector<sgNodeID_t> nodes);
    /**
     * Adds non-anchor intermediate nodes to a line by finding the transitive path between anchors in a multi_ldg.
     * If multiple paths exist, use evidence to decide or avoid filling.
     * @param line
     * @param multi_ldg
     * @param min_links
     * @return
     */
    //std::vector<sgNodeID_t> add_intermediate_nodes(std::vector<sgNodeID_t> line, const LinkageDiGraph & multi_ldg, int min_links=3);

    WorkSpace &ws;
    //std::vector<LinkageDiGraph> multi_linkage_graphs;
    //std::vector<LinkageDiGraph> anchor_linkage_graphs;

    std::vector<bool> selected_nodes;
    std::vector<bool> frontier_nodes;
};


#endif //BSG_LINKAGEUNTANGLER_HPP
