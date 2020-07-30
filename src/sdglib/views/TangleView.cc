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