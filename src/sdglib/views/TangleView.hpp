//
// Created by Bernardo Clavijo (EI) on 30/07/2020.
//

#pragma once

#include <sdglib/views/NodeView.hpp>

/**
 * A tangle is an unresolved component of the graph.
 * Internal nodes are only connected to other internal nodes and to frontiers, where the frontiers are directional outwards.
 * All connections from/to a frontier (in directional manner) are included in the Tangle, even other Frontiers that are only connected to frontiers.
 *
 * is_complete checks the sanity of the tangle connectivity, without any knowledge of why a frontier is a frontier.
 */
class TangleView {
public:
    TangleView(DistanceGraph const * dg, const std::vector<sgNodeID_t> & _frontiers,const std::vector<sgNodeID_t> & _internals);
    bool is_complete();
    void dump_subgraph(int radious=2);

    std::string classify_tangle() const;

    std::vector<NodeView> frontiers;
    std::vector<NodeView> internals;
    const DistanceGraph * dg;
};
