//
// Created by Bernardo Clavijo (EI) on 2019-06-26.
//

#include "NodeView.hpp"
#include <sdglib/graph/SequenceDistanceGraph.hpp>

DistanceGraph NodeView::graph() const {
    return *dg;
}

const std::string NodeView::sequence() const {
    return dg->sdg.get_node_sequence(id);
}

const uint64_t NodeView::size() const {
    return dg->sdg.get_node_size(id);
}

const std::vector<NodeDistanceView> NodeView::next() const {
    auto fwl=dg->get_fw_links(id);
    std::vector<NodeDistanceView> r;
    r.reserve(fwl.size());
    for (auto &l:fwl) {
        r.emplace_back(dg->get_nodeview(l.dest),l.dist,l.support);
    }
    return r;
}

const std::vector<NodeDistanceView> NodeView::prev() const {
    auto bwl=dg->get_bw_links(id);
    std::vector<NodeDistanceView> r;
    r.reserve(bwl.size());
    for (auto &l:bwl) {
        r.emplace_back(dg->get_nodeview(-l.dest),l.dist,l.support);
    }
    return r;
}

std::ostream &operator<<(std::ostream &os, const NodeView &nv) {
    os << "< NodeView: "<<nv.id<<" in "<<nv.dg->name<<" >";
    return os;
}

std::ostream &operator<<(std::ostream &os, const NodeDistanceView &ndv) {
    os << "< NodeDistanceView: "<<ndv.dist<<"bp to "<<ndv.node_view.node_id()<<" >";
    return os;
}