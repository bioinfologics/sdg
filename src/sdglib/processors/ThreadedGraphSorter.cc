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

bool ThreadedGraphSorter::pop_node(sgNodeID_t node_id, uint64_t read){
    auto nv = dg.get_nodeview(node_id);
    std::vector<LinkView> ln;
    std::vector<LinkView> lp;
    // Get reads supporting node connections
    for (const auto& lv: nv.next()){
        if (lv.support().id == read) ln.push_back(lv);
    }
    for (const auto& lv: nv.prev()){
        if (lv.support().id == read) lp.push_back(lv);
    }

    // Tip case
    if (ln.empty() and lp.size()==1){
        dg.remove_link(-lp[0].node().node_id(), node_id, lp[0].distance(), lp[0].support());
        return true;
    } else if (ln.size() == 1 and lp.empty()) {
        dg.remove_link(-node_id, lp[0].node().node_id(), lp[0].distance(), lp[0].support());
        return true;
    }

    // TODO ask here!! --> node is tip or it's connected to itself???
    if (ln.size()!=1 or lp.size()!=1 or llabs(ln[0].node().node_id()) == llabs(node_id) or llabs(lp[0].node().node_id()) == llabs(node_id))
        return false;

    // Pops the node from the line
    uint32_t td = ln[0].distance()+nv.size()+lp[0].distance();
    dg.add_link(-lp[0].node().node_id(), ln[0].node().node_id(), td, lp[0].support());
    dg.remove_link(-node_id,ln[0].node().node_id(),ln[0].distance(),ln[0].support());
    dg.remove_link(-lp[0].node().node_id(),node_id,lp[0].distance(),lp[0].support());
    return true;
}

