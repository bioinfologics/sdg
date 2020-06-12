//
// Created by Bernardo Clavijo (EI) on 20/05/2020.
//

#include <sdglib/views/NodeView.hpp>
#include "GraphPatcher.hpp"

uint64_t GraphPatcher::find_tips_to_reconnect(int min_paths) {
    for (auto fnv:ws.sdg.get_all_nodeviews())
        for (auto nv: {fnv,fnv.rc()}) {
            if (nv.next().empty()) {
                std::unordered_map<sgNodeID_t,uint64_t> votes;
                for (auto &fwp: ws.paired_reads_datastores[0].mapper.all_paths_fw(nv.node_id())) {
                    if (fwp[0]==0) ++votes[fwp[1]];
                    else ++votes[fwp[0]];
                }
                std::vector<sgNodeID_t> rg;
                rg.emplace_back(nv.node_id());
                for (auto &v:votes) {
                    if (v.second>min_paths) rg.emplace_back(-v.first);
                }
                if (rg.size()>1) {
                    reconnection_groups.emplace_back(rg);
                }
            }
        }
    return reconnection_groups.size();

}