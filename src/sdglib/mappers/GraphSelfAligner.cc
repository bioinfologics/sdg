//
// Created by Bernardo Clavijo (EI) on 2019-08-22.
//

#include <sdglib/indexers/NKmerIndex.hpp>
#include <sdglib/mappers/SequenceMapper.hpp>
#include "GraphSelfAligner.hpp"
#include <sdglib/views/NodeView.hpp>

void GraphSelfAligner::self_align() {
    matches.clear();
    matches.resize(dg.sdg.nodes.size());
    SequenceMapper sm(dg,15);
    sm.update_graph_index(6);
    for (auto &nv:dg.get_all_nodeviews()){
        for (auto &m:sm.map_sequence(nv.sequence().c_str(),nv.node_id())){
            if (llabs(m.node)!=nv.node_id()){
                std::cout<<"Adding match from "<<nv<<" to matches["<<llabs(nv.node_id())<<"] :"<<m<<std::endl;
                matches[llabs(nv.node_id())].emplace_back(nv.node_id(),m.node,m.qStart,m.qEnd,m.nStart,m.nEnd,m.score);
            }
        }

    }
    std::sort(matches.begin(),matches.end());
}