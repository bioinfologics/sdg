//
// Created by Gonzalo Garcia (EI) on 2020-11-10.
//

#include "ThreadedGraphSorter.h"

std::vector<uint64_t > ThreadedGraphSorter::rids_from_node(NodeView nv){
    // TODO: report that the rids are uint64_t but the ids in support are int64_t
    std::unordered_set<uint64_t > rids;
    for (const auto& lv: nv.next()){
        rids.insert((uint64_t)lv.support().id);
    }
    for (const auto& lv: nv.prev()){
        rids.insert((uint64_t)lv.support().id);
    }
    return std::vector<uint64_t >(rids.begin(), rids.end());
}

uint64_t ThreadedGraphSorter::shared_reads(NodeView nv1, NodeView nv2){
    uint64_t shared=0;
    auto r1 = rids_from_node(nv1);
    auto r2 = rids_from_node(nv2);
    std::sort(r1.begin(), r1.end());
    std::sort(r2.begin(), r2.end());
    
    int pr1=0;
    int pr2=0;
    while (pr1<r1.size() and pr2<r2.size()){
        if (r1[pr1]<r2[pr2]){
            pr1++;
        } else if (r1[pr1]>r2[pr2]){
            pr2++;
        } else {
            shared++;
            pr1++;
            pr2++;
        }
    }
    return shared;
}