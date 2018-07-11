//
// Created by Bernardo Clavijo (EI) on 28/05/2018.
//

#ifndef BSG_LINKAGEUNTANGLER_HPP
#define BSG_LINKAGEUNTANGLER_HPP


#include <sglib/WorkSpace.hpp>
#include <sglib/LinkageDiGraph.hpp>

class LinkageUntangler {
public:

    explicit LinkageUntangler(WorkSpace & _ws): ws(_ws) { clear_node_selection();};

    //Node selection methods
    void clear_node_selection();
    void report_node_selection();
    void select_nodes_by_size_and_ci( uint64_t min_size, float min_ci, float max_ci);
    std::set<std::pair<sgNodeID_t, sgNodeID_t >> get_HSPNPs(uint64_t min_size, float min_ci, float max_ci);
    void select_nodes_by_HSPNPs(uint64_t min_size, float min_ci, float max_ci);
    void select_frontiers_by_size_and_ci();

    //Linkage creation methods (work on selected nodes)
    LinkageDiGraph make_topology_linkage(int radius);
    LinkageDiGraph make_paired_linkage(int min_reads);
    LinkageDiGraph make_tag_linkage(int min_tags,float end_perc=.3);

    //Linkage filtering methods
    LinkageDiGraph filter_linkage_to_hspnp_duos( uint64_t min_size, float min_ci, float max_ci, const LinkageDiGraph & ldg);

    //Graph untangling methods
    void expand_trivial_repeats(const LinkageDiGraph &);
    void expand_linear_regions(const LinkageDiGraph &);
    void expand_linear_regions_skating(const LinkageDiGraph &, int max_lines=0);
    void linear_regions_tag_local_assembly(const LinkageDiGraph & ldg, uint8_t k, int min_cvg, int max_lines, uint64_t min_nodes, uint64_t min_total_size, bool count_tag_cvg=false);
    //Problem localisation methods


    WorkSpace &ws;
    std::vector<bool> selected_nodes;
    std::vector<bool> frontier_nodes;
};


#endif //BSG_LINKAGEUNTANGLER_HPP
