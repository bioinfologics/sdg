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
 * Can select nodes which are to be used for linkage.
 */
class LinkageUntangler {
    WorkSpace &ws;
public:

    explicit LinkageUntangler(WorkSpace & _ws): ws(_ws) { clear_node_selection();};

    //Node selection methods
    void clear_node_selection();
    void report_node_selection();
    void select_nodes_by_size_and_ci( uint64_t min_size, float min_ci, float max_ci);
    std::set<std::pair<sgNodeID_t, sgNodeID_t >> get_HSPNPs(uint64_t min_size, float min_ci, float max_ci);
    void select_nodes_by_HSPNPs(uint64_t min_size, float min_ci, float max_ci);
    //void select_frontiers_by_size_and_ci();

    //Linkage creation methods (work on selected nodes)
    LinkageDiGraph make_topology_linkage(int radius);
    LinkageDiGraph make_paired_linkage(int min_reads);
    LinkageDiGraph make_tag_linkage(int min_tags,float end_perc=.3);
    LinkageDiGraph make_longRead_linkage();

    //Linkage filtering methods
    LinkageDiGraph filter_linkage_to_hspnp_duos( uint64_t min_size, float min_ci, float max_ci, const LinkageDiGraph & ldg);

    //Graph untangling methods
    void expand_trivial_repeats(const LinkageDiGraph &);
    void expand_linear_regions(const LinkageDiGraph &);
    void expand_linear_regions_skating(const LinkageDiGraph &, int max_lines=0);

    //Problem localisation methods


    std::vector<bool> selected_nodes;
    std::vector<bool> frontier_nodes;
};


#endif //BSG_LINKAGEUNTANGLER_HPP
