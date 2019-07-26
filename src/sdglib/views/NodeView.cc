//
// Created by Bernardo Clavijo (EI) on 2019-06-26.
//

#include "NodeView.hpp"
#include <sdglib/graph/SequenceDistanceGraph.hpp>
#include <sdglib/workspace/WorkSpace.hpp>

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
    std::sort(r.begin(),r.end());
    return r;
}

const std::vector<NodeDistanceView> NodeView::prev() const {
    auto bwl=dg->get_bw_links(id);
    std::vector<NodeDistanceView> r;
    r.reserve(bwl.size());
    for (auto &l:bwl) {
        r.emplace_back(dg->get_nodeview(-l.dest),l.dist,l.support);
    }
    std::sort(r.begin(),r.end());
    return r;
}

std::vector<uint16_t> NodeView::kmer_coverage(std::string kcovds_name, std::string kcovds_count_name) const {
    return dg->sdg.ws.get_kmer_counter(kcovds_name).project_count(kcovds_count_name,sequence());
}

std::vector<uint16_t> NodeView::kmer_coverage(int kcovds_idx, int kcovds_count_idx) const {
    return dg->sdg.ws.kmer_counters[kcovds_count_idx].project_count(kcovds_count_idx,sequence());
}

std::ostream &operator<<(std::ostream &os, const NodeView &nv) {
    os << "< NodeView: "<<nv.id<<" in "<<nv.dg->name<<" >";
    return os;
}

std::ostream &operator<<(std::ostream &os, const NodeDistanceView &ndv) {
    os << "< NodeDistanceView: "<<ndv.dist<<"bp to "<<ndv.node_view.node_id()<<" >";
    return os;
}