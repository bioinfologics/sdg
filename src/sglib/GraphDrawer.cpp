//
// Created by Ben Ward (EI) on 08/02/2018.
//

#include "GraphDrawer.h"

GraphDrawer::GraphDrawer(const SequenceGraph& _sg) : sg(_sg) {
    // Create canonical set of links.
    for (auto& links : sg.links) {
        for (auto& l : links) {
            canonical_edge_pool.emplace(l.make_canonical());
        }
    }
};

std::string GraphDrawer::render_node_as_dot_string(const unsigned long i, const Node& n) const {
    std::stringstream dotstring;
    dotstring << "seqnode" << i << " [label=\"<source> + | Node " << i;
    dotstring << " | " << n.sequence.size() << " bp |<sink> - \"];";
    return dotstring.str();
}

std::string GraphDrawer::make_dot_edge(const sgNodeID_t from, const sgNodeID_t to) const {
    std::stringstream dotstring;
    std::string insign = from > 0 ? "source" : "sink";
    std::string outsign = to > 0 ? "source" : "sink";
    dotstring << "seqnode" << std::abs(from) << ':' << insign;
    dotstring << " -> ";
    dotstring << "seqnode" << std::abs(to) << ':' << outsign;
    return dotstring.str();
}

std::string GraphDrawer::make_dot_edge(const Link& l) const {
    return make_dot_edge(l.source, l.dest);
}

std::string GraphDrawer::render_edge_as_dot_string(const Link& l) const {
    return make_dot_edge(l) + ";";
}

void GraphDrawer::write_dot_header(std::ofstream& gvfile) const {
    gvfile << "digraph seqgraph {" << std::endl;
    gvfile << "\t// Apply the defaults" << std::endl;
    gvfile << "\tnode [shape=Mrecord];" << std::endl;
}

void GraphDrawer::write_nodes_to_dot(std::ofstream& gvfile) const {
    gvfile << "\t// Create all nodes with default values" << std::endl;
    for (auto& n : sg.nodes) {
        auto nodestr = render_node_as_dot_string(&n - &sg.nodes[0], n);
        gvfile << '\t' << nodestr << std::endl;
    }
}

void GraphDrawer::write_non_path_links_to_dot(std::ofstream& gvfile) const {
    gvfile << "\t// Create all links with default values" << std::endl;
    gvfile << "\tsubgraph allconnections {" << std::endl;
    gvfile << "\t\tedge [dir=none, style=dotted];" << std::endl;
    for (auto& l : canonical_edge_pool) {
        auto linkstr = render_edge_as_dot_string(l);
        gvfile << "\t\t" << linkstr << std::endl;
    }
    gvfile << "\t}" << std::endl;
}

void GraphDrawer::render_graph_as_dot(std::ofstream& gvfile) const {
    write_dot_header(gvfile);
    write_nodes_to_dot(gvfile);
    write_non_path_links_to_dot(gvfile);
    gvfile << "\t// Adding paths as subgraphs" << std::endl;
    for (auto& path : paths) {
        gvfile << path << std::endl;
    }
    gvfile << "\t Adding extra attributes" << std::endl;
    for (const auto& attribute : extra_attributes) {
        gvfile << attribute << std::endl;
    }
    gvfile << '}' << std::endl;
    gvfile.close();
}

void GraphDrawer::style_nodes(const std::vector<sgNodeID_t>& nodes, const std::string& colouring) {
    std::stringstream nodestring;
    nodestring << '\t';
    auto node = nodes.begin();
    for (; node != std::prev(nodes.end()); ++node) {
        nodestring << "seqnode" << std::abs(*node) << ", ";
    }
    nodestring << "seqnode" << std::abs(*node) << " [color=" << colouring << "];";
    extra_attributes.emplace_back(nodestring.str());
}

void GraphDrawer::add_path(const std::string& name, const SequenceGraphPath& path, const std::string& colouring, const bool colournodes) {
    if (colournodes) {
        style_nodes(path.nodes, colouring);
    }
    auto edges = path.collect_links();
    std::stringstream edgestring;
    if(edges.size() != 0) {
        edgestring << "\tsubgraph " << name << "{" << std::endl << "\t\tedge [color=" << colouring << "]" << std::endl;
        for (const auto& edge : edges) {
            auto canonical_edge = edge.make_canonical();
            auto edge_location = canonical_edge_pool.find(canonical_edge);
            if(edge_location == canonical_edge_pool.end()) {
                std::cout << "ERROR: ADDING PATH USING EDGE NOT IN THE CANONICAL EDGE POOL" << std::endl;
            } else {
                canonical_edge_pool.erase(edge_location);
            }
            edgestring << "\t\t" << make_dot_edge(edge) << std::endl;
        }
        edgestring << "\t}" << std::endl;
        paths.emplace_back(edgestring.str());
    }
}

std::vector<std::string> GraphDrawer::brew_colours(uint64_t n) const {
    std::vector<std::string> hsv_vals;
    auto angle_per_slice = 360.0 / n;
    for (auto i = 0; i < n; i++) {
        auto hsv = std::to_string(i * angle_per_slice) + "1.0 1.0";
        hsv_vals.emplace_back(hsv);
    }
    return hsv_vals;
}