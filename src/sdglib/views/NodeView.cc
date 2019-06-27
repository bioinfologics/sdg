//
// Created by Bernardo Clavijo (EI) on 2019-06-26.
//

#include "NodeView.hpp"
#include <sdglib/graph/SequenceDistanceGraph.hpp>

std::string NodeView::sequence() const {
    return dg.sdg.get_node_sequence(node_id);
}

uint64_t NodeView::size() const {
    return dg.sdg.get_node_size(node_id);
}

std::vector<NodeDistanceView> NodeView::next() const {
    auto fwl=dg.get_fw_links(node_id);
    std::vector<NodeDistanceView> r;
    r.reserve(fwl.size());
    for (auto &l:fwl) {
        r.emplace_back(NodeView(dg,l.dest),l.dist,l.support);
    }
    return r;
}

std::vector<NodeDistanceView> NodeView::prev() const {
    auto bwl=dg.get_bw_links(node_id);
    std::vector<NodeDistanceView> r;
    r.reserve(bwl.size());
    for (auto &l:bwl) {
        r.emplace_back(NodeView(dg,-l.dest),l.dist,l.support);
    }
    return r;
}

std::ostream &operator<<(std::ostream &os, const NodeView &nv) {
    os << "< NodeView: "<<nv.node_id<<" in "<<nv.dg.name<<" >";
    return os;
}

std::ostream &operator<<(std::ostream &os, const NodeDistanceView &ndv) {
    os << "< NodeDistanceView: "<<ndv.distance<<"bp to "<<ndv.node.node_id<<" >";
    return os;
}

NodeView& NodeView::operator=(NodeView const &o) {
    if (&o == this) return *this;
    dg=o.dg;
    node_id=o.node_id;
}

NodeDistanceView& NodeDistanceView::operator=(NodeDistanceView const &o) {
    if (&o == this) return *this;
    node=o.node;
    distance=o.distance;
    support=o.support;
}