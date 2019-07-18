//
// Created by Bernardo Clavijo (EI) on 28/05/2018.
//

#pragma once


#include <sdglib/workspace/WorkSpace.hpp>
#include <sdglib/graph/DistanceGraph.hpp>

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

    explicit LinkageUntangler(const DistanceGraph & _dg): dg(_dg) { deselect_all();};

    explicit LinkageUntangler(const SequenceDistanceGraph &_dg):LinkageUntangler(static_cast<const DistanceGraph &>(_dg)){};

    //==== Methods currently in use =====

    //Node selection methods
    void deselect_all();
    void select_all();
    void report_selection();
    void select_by_size(uint64_t min_size,uint64_t max_size=0);

    //std::set<std::pair<sgNodeID_t, sgNodeID_t >> get_HSPNPs(uint64_t min_size, float min_ci, float max_ci);

    /**
    * Anchors are chosen by having all their strongly connected neighbours (i.e. those with min_links connections)
    * coherently connected among them in order with at least min_transitive_links.
    * @param min_links
    * @param min_transitive_links
    */
    void select_linear_anchors(int min_links=3, int min_transitive_links=2);

    //Multi-Linkage creation methods: multiple evidence-supported links between any nodes



    //Linkage creation methods: aggregated links between selected nodes (anchors)

    DistanceGraph make_nextselected_linkage(int min_links=3);

    //DistanceGraph make_and_simplify_anchor_linkage(const DistanceGraph &multi_ldg, int anchor_link=5, int transitive_links=3);



    //==== DEPRECATED methods, do NOT use in new projects =====


    //void select_nodes_by_HSPNPs(uint64_t min_size, float min_ci, float max_ci);



    //Linkage improving/filtering methods
    //DistanceGraph filter_linkage_to_hspnp_duos( uint64_t min_size, float min_ci, float max_ci, const DistanceGraph & ldg);


    //Graph untangling/modification/local assembly methods
//    void expand_trivial_repeats(const DistanceGraph &);
//    void expand_linear_regions(const DistanceGraph &);
//    void expand_linear_regions_skating(const DistanceGraph &, int max_lines=0);
//    void linear_regions_tag_local_assembly(const DistanceGraph & ldg, uint8_t k, int min_cvg, int max_lines, uint64_t min_nodes, uint64_t min_total_size, bool count_tag_cvg=false);
//
//    //Problem localisation methods
//    void fill_linkage_line(std::vector<sgNodeID_t> nodes);
    /**
     * Adds non-anchor intermediate nodes to a line by finding the transitive path between anchors in a multi_ldg.
     * If multiple paths exist, use evidence to decide or avoid filling.
     * @param line
     * @param multi_ldg
     * @param min_links
     * @return
     */
    //std::vector<sgNodeID_t> add_intermediate_nodes(std::vector<sgNodeID_t> line, const DistanceGraph & multi_ldg, int min_links=3);

    const DistanceGraph & dg;

    std::vector<bool> selected_nodes;
};

