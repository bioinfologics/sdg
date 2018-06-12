//
// Created by Bernardo Clavijo (EI) on 25/05/2018.
//

#ifndef BSG_LINKAGEDIGRAPH_HPP
#define BSG_LINKAGEDIGRAPH_HPP


#include "SequenceGraph.hpp"

class LinkageDiGraph {
public:
    LinkageDiGraph(SequenceGraph & _sg): sg(_sg){};

    void add_link( sgNodeID_t source, sgNodeID_t dest, int32_t d);
    void add_links(const LinkageDiGraph &other);

    void remove_link(sgNodeID_t source, sgNodeID_t dest);

    std::vector<Link> get_fw_links( sgNodeID_t n) const;
    std::vector<Link> get_bw_links( sgNodeID_t n) const;
    std::set<sgNodeID_t> fw_reached_nodes(sgNodeID_t n, int radius) const;

    void remove_transitive_links(int radius);
    void report_connectivity();
    void solve();


    SequenceGraph & sg;
    std::vector<std::vector<Link>> links;

};
#endif //BSG_LINKAGEDIGRAPH_HPP
