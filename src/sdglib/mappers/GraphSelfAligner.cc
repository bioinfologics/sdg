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
    SequenceMapper sm(dg,31,6);
    auto nvs=dg.get_all_nodeviews();
#pragma omp parallel for  schedule(static,200)
    for (auto i=0;i<nvs.size();++i){
        auto &nv=nvs[i];
        for (auto &m:sm.map_sequence(nv.sequence().c_str(),nv.node_id())){
            if (llabs(m.node)!=nv.node_id()){
                matches[llabs(nv.node_id())].emplace_back(nv.node_id(),m.node,m.qStart,m.qEnd,m.nStart,m.nEnd,m.score);
            }
        }
        std::sort(matches[llabs(nv.node_id())].begin(),matches[llabs(nv.node_id())].end());
    }
}