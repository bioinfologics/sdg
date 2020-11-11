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

bool ThreadedGraphSorter::multipop_node(sgNodeID_t node_id, uint64_t read){
    // TODO: what for no lo entiendo??
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

    if (ln.size()!=1 or lp.size()!=1 or llabs(ln[0].node().node_id())==llabs(node_id) or llabs(lp[0].node().node_id())==llabs(node_id)) {
        for (const auto& p: lp){
            for (const auto& n: ln){
                std::cout<<p.node() <<" "<< n.node()<< " " << shared_reads(p.node(), n.node())<<std::endl;
            }
        }
        return false;
    }
    return pop_node(node_id, read);
}

void ThreadedGraphSorter::write_connected_nodes_graph(std::string filename){
    std::vector<sgNodeID_t > sn;

    // TODO: for some reason the include_diconnected parameters is giving an error!
//    for (const auto& nv: dg.get_all_nodeviews(include_disconnected=false)){
    for (const auto& nv: dg.get_all_nodeviews()){
        if (!dg.links[llabs(nv.node_id())].empty()){
            sn.push_back(nv.node_id());
        }
    }
    dg.write_to_gfa1(filename, sn);
}

std::pair<std::map<sgNodeID_t , int64_t >, std::vector<sgNodeID_t >> ThreadedGraphSorter::sort_cc(std::vector<sgNodeID_t> cc){
    //uses relative position propagation to create a total order for a connected component
    //    returns a dict of nodes to starting positions, and a sorted list of nodes with their positions

    // Next nodes to process
    std::vector<sgNodeID_t > next_nodes;
    // Map with the nodes positions
    std::map<sgNodeID_t , int64_t > node_pos;

    // check no node on the CC appears in both directions.
    for (const auto&x: cc){
        if (std::find(cc.begin(), cc.end(), -x)!=cc.end()){
            return std::make_pair(node_pos, next_nodes);
        }
    }

    // now start from each node without prev as 0
    for (const auto& x: cc){
        if (dg.get_nodeview(x).prev().empty()){
            next_nodes.push_back(x);
        }
        // Nodes with no prev are the inital nodes in the partial order
        node_pos[x]=0;
    }

    // Propagate position FW
    std::vector<std::vector<sgNodeID_t >> discovered_paths;
    while (!next_nodes.empty()) {
        std::vector<sgNodeID_t > new_next_nodes;
        for (const auto& nid: next_nodes){
            for (const auto& l: dg.get_nodeview(nid).next()){
                auto nid_node_position = node_pos[nid];
                auto nid_node_size = dg.get_nodeview(nid).size();
                auto linked_node_id = l.node().node_id();
                auto linked_node_position = node_pos[linked_node_id];
                // TODO: does this work like this --> if nid node ends after the linked node the next node is moved fw and the linked node is added to the list to place
                if (nid_node_position+nid_node_size+l.distance()>linked_node_position){
                    node_pos[linked_node_id] = nid_node_position+nid_node_size+l.distance();
                    new_next_nodes.push_back(linked_node_id);
                }
            }
        }
        discovered_paths.push_back(next_nodes);
        next_nodes=new_next_nodes;
    }
    // move all nodes with no next to the end by asigning the max position to all nodes with no next
    int64_t lastp=0;
    for (const auto& mi: node_pos){
        if (mi.second > lastp) lastp=mi.second;
    }
    for (const auto& nid: cc){
        if (dg.get_nodeview(nid).next().empty()){
            node_pos[nid]=lastp;
        }
    }

    // move every node up to its limit fw
    bool updated=true;
    while (updated){
        updated=false;
        //sort node_pos by the second element
        // TODO: This sorted_node_pos could be copies once to the vector outside the loop and then work with the vector until return
        std::vector<std::pair<sgNodeID_t , int64_t>> sorted_node_pos;
        std::copy(node_pos.begin(), node_pos.end(), std::back_inserter< std::vector<std::pair<sgNodeID_t , int64_t>>>(sorted_node_pos));
        std::sort(sorted_node_pos.begin(), sorted_node_pos.end(), [](std::pair<sgNodeID_t , int64_t> elem1 ,std::pair<sgNodeID_t , int64_t> elem2){return elem1.second < elem2.second;});

        for (const auto& mp: sorted_node_pos){
            auto nid = mp.first;
            auto p = mp.second;
            auto nid_node_size = dg.get_nodeview(nid).size();

            if (!dg.get_nodeview(nid).next().empty()){
                // get min distance between the nid and the others, means the next node
                int64_t np = std::numeric_limits<int64_t >::max;
                for (const auto& l: dg.get_nodeview(nid).next()){
                    auto linked_node_id = l.node().node_id();
                    if (node_pos[linked_node_id]-l.distance()-nid_node_size < np)
                        np = node_pos[linked_node_id]-l.distance()-nid_node_size;
                }
                if (np!=node_pos[nid]){
                    node_pos[nid]=np;
                    updated=true;
                }

            }
        }
    }
    // move every node up to its limit bw
    updated=true;
    while(updated){
        updated=false;
        // TODO: This sorted_node_pos could be copies once to the vector outside the loop and then work with the vector until return
        std::vector<std::pair<sgNodeID_t , int64_t>> sorted_node_pos;
        std::copy(node_pos.begin(), node_pos.end(), std::back_inserter< std::vector<std::pair<sgNodeID_t , int64_t>>>(sorted_node_pos));
        std::sort(sorted_node_pos.begin(), sorted_node_pos.end(), [](std::pair<sgNodeID_t , int64_t> elem1 ,std::pair<sgNodeID_t , int64_t> elem2){return elem1.second < elem2.second;});

        for (const auto& mp: sorted_node_pos){
            auto nid = mp.first;
            auto p = mp.second;
            auto nid_node_size = dg.get_nodeview(nid).size();

            if (!dg.get_nodeview(nid).next().empty()){
                // get min distance between the nid and the others, means the next node
                int64_t np = 0;
                for (const auto& l: dg.get_nodeview(nid).next()){
                    auto linked_node_id = l.node().node_id();
                    if (node_pos[linked_node_id]+l.distance()+nid_node_size > np)
                        np = node_pos[linked_node_id]+l.distance()+nid_node_size;
                }
                if (np!=node_pos[nid]){
                    node_pos[nid]=np;
                    updated=true;
                }
            }
        }
    }
    return std::make_pair(node_pos, next_nodes);
}

