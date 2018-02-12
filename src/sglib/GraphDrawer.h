//
// Created by Ben Ward (EI) on 08/02/2018.
//

#ifndef BSG_GRAPHDRAWER_H
#define BSG_GRAPHDRAWER_H

#include "SequenceGraph.hpp"
#include <fstream>
#include <sstream>
#include <set>

class GraphDrawer {
private:
    const SequenceGraph& sg;
    std::set<Link> canonical_edge_pool;
    std::vector<std::string> extra_attributes;
    std::vector<std::string> paths;
    std::vector<SequenceGraphPath> sgpaths;

    std::string render_node_as_dot_string(unsigned long i, const Node& n) const;
    std::string make_dot_edge(sgNodeID_t from, sgNodeID_t to) const;
    std::string make_dot_edge(const Link& l) const;
    std::string render_edge_as_dot_string(const Link& l) const;

    void write_dot_header(std::ofstream& gvfile) const;
    void write_nodes_to_dot(std::ofstream& gvfile) const;
    void write_non_path_links_to_dot(std::ofstream& gvfile) const;

public:
    explicit GraphDrawer(const SequenceGraph& _sg);

    std::vector<std::string> brew_colours(uint64_t n) const;
    void render_graph_as_dot(std::ofstream& file) const;
    void add_path(const std::string& name, const SequenceGraphPath& path, const std::string& colouring, bool colournodes);
    void style_nodes(const std::vector<sgNodeID_t>& nodes, const std::string& colouring);
};











#endif //BSG_GRAPHDRAWER_H
