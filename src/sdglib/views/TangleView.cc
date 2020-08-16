//
// Created by Bernardo Clavijo (EI) on 30/07/2020.
//

#include "TangleView.hpp"

TangleView::TangleView(DistanceGraph const * _dg, const std::vector<sgNodeID_t> &_frontiers, const std::vector<sgNodeID_t> &_internals) : dg(_dg){
    frontiers.reserve(_frontiers.size());
    for (auto &f:_frontiers) frontiers.emplace_back(NodeView(dg,f));
    internals.reserve(_internals.size());
    for (auto &i:_internals) internals.emplace_back(NodeView(dg,i));
}

bool TangleView::is_complete() {
    //set of all nodes at play (frontiers are only included in their outgoing direction)

    //check all nexts/prevs for internals are in the set

    //check all prevs for frontiers are in the set
    return false;
}


void TangleView::dump_subgraph(int radious){

}

std::string TangleView::classify_tangle() const {

    if (frontiers.size() == 0) {
        return "debris";
    }
    if (frontiers.size() == 1) {

        std::unordered_set<sgNodeID_t> nids;
        for (const auto &nv: internals) nids.emplace(nv.node_id());
        std::unordered_set<sgNodeID_t> seen_nids;

        std::vector<NodeView> nexts{frontiers[0].rc()};
        bool abort = false;
        while (nexts.size() > 0 and !abort) {
            std::vector<NodeView> new_nexts;
            for (const auto &nv: nexts) {
                for (const auto &nnl: nv.next()) {
                    if (std::find(nids.begin(), nids.end(), std::abs(nnl.node().node_id())) != nids.end()) {
                        abort = true;
                        break;
                    }
                    if (std::find(seen_nids.begin(), seen_nids.end(), std::abs(nnl.node().node_id())) !=
                        seen_nids.end()) {
                        seen_nids.emplace(nnl.node().node_id());
                        new_nexts.push_back(nnl.node());
                    }
                }
                if (abort) break;
            }
            nexts = new_nexts;
        }
        if (!abort) return "complex tip";
    }
    if (internals.size() == 1) {
        auto nv = internals[0];
        if (frontiers.size() == 2 and nv.prev().size() + nv.next().size() == 1) return "tip";

        std::unordered_set<sgNodeID_t> involved_nodes_set{std::abs(nv.node_id())};
        for (const auto &tnv: nv.prev()) {
            involved_nodes_set.emplace(std::abs(tnv.node().node_id()));
        }
        for (const auto &tnv: nv.next()) {
            involved_nodes_set.emplace(std::abs(tnv.node().node_id()));
        }

        if (frontiers.size() == 2 * nv.prev().size() and nv.prev().size() == nv.next().size() and
            involved_nodes_set.size() == nv.prev().size() * 2 + 1) {
            return "repeat " + std::to_string(nv.prev().size()) + ":" + std::to_string(nv.next().size());
        }
    }
    if (internals.size() == 2 and frontiers.size() == 2) {
        auto i0p = internals[0].parallels();
        if (i0p.size() == 1 and i0p[0] == internals[1]) {
            return "bubble";
        } else {
            auto i0n = internals[0].next();
            auto i0p = internals[0].prev();
            auto i1n = internals[1].next();
            auto i1p = internals[1].prev();

            if (i0n.size() == 2 and i0p.size() == 2 and i1n.size() == 1 and i1p.size() == 1 and
                i1n[0].node() == i1p[0].node()) {
                return "loop";
            }
            if (i1n.size() == 2 and i1p.size() == 2 and i0n.size() == 1 and i0p.size() == 1 and
                i0n[0].node() == i0p[0].node()) {
                return "loop";
            }
        }
    }
    return "unclassified";
};