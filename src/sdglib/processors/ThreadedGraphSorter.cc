//
// Created by Gonzalo Garcia (EI) on 2020-11-10.
//

#include "ThreadedGraphSorter.h"
#include <sdglib/graph/ReadThreadsGraph.hpp>
#include <atomic>
enum Happiness {unknown=-1,unhappy=0,happy=1};

std::array<uint64_t,3> assess_node_happiness(sgNodeID_t nid, const std::unordered_map<sgNodeID_t , uint32_t> & order, const DistanceGraph& trg_nt){

    // Check that nid or -nid is not in order
    if (order.find(nid)==order.end())
        nid*=-1;
    if (order.find(nid)==order.end())
        return {0,0,0};

    auto nv=trg_nt.get_nodeview(nid);
    std::unordered_set<uint64_t > rids;
    for (const auto &lv: nv.prev())
        rids.insert(lv.support().id);
    for (const auto &lv: nv.next())
        rids.insert(lv.support().id);

    auto npos = order[nid];
    uint64_t disconnected=0;
    uint64_t happy=0;
    uint64_t unhappy=0;
    for (const auto &rid: rids){
        Happiness happy_fw=Happiness::unknown;
        std::vector<LinkView > nlvs; //TODO: replace for a single nodeview with nid=0 if not yet filled.
        for (const auto& lv: nv.next()){
            if (lv.support().id==rid){
                nlvs.push_back(lv);
            }
        }
        while (!nlvs.empty()){
            auto nnv=nlvs[0].node();
            if (order.find(nnv.node_id())!=order.end()){
                if (order[nnv.node_id()]>npos){
                    happy_fw=Happiness::happy;
                } else {
                    happy_fw=Happiness::unhappy;
                    break;
                }
            }
            nlvs.clear();
            for (const auto& l: nnv.next()){
                if (l.support().id == rid){
                    nlvs.push_back(l);
                }
            }
        }
        int happy_bw=Happiness::unknown;
        std::vector<LinkView > pnvs;
        for (const auto& l: nv.prev()){
            if (l.support().id==rid){
                pnvs.push_back(l);
            }
        }
        while (!pnvs.empty()){
            auto pnv = pnvs[0].node();
            if (order.find(pnv.node_id())!=order.end()){
                if (order[pnv.node_id()]<npos){
                    happy_bw=Happiness::happy;
                } else {
                    happy_bw=Happiness::unhappy;
                    break;
                }
            }
            pnvs.clear();
            for (const auto& l: pnv.prev()){
                if (l.support().id == rid){
                    pnvs.push_back(l);
                }
            }
        }

        if (happy_fw==Happiness::unhappy or happy_bw==Happiness::unhappy){
            //At least one end is unhappy, glass half empty !
            unhappy++;
        }
        else if (happy_fw==Happiness::unknown and happy_bw==Happiness::unknown){
            //We don't really know how we feel
            disconnected++;
        }
        else {
            //At least one end is happy, and no ends are unhappy, glass half full !
            happy++;
        }
    }
    return {happy, unhappy, disconnected};
}

int pop_node_from_all(DistanceGraph &dg, sgNodeID_t node_id){
    auto nv = dg.get_nodeview(node_id);
    std::vector<uint64_t > rids;
    for (const auto& rid: nv.next())
        rids.push_back(rid.support().id);
    for (const auto& rid: nv.prev())
        rids.push_back(rid.support().id);

    int popped_count=0;
    for (const auto &rid: rids){
        if (pop_node(dg, node_id, rid))
            popped_count++;
    }
    return popped_count;
}

std::vector<NodePosition> make_thread_happy(const std::vector<NodePosition> &thread,const DistanceGraph & trg, int max_unhappy, float disconnection_rate){

    std::unordered_map<int64_t,std::vector<std::pair<int64_t,uint16_t>>> nodepos_in_reads;//for every read present in any node, get a list of all nodes in it
    for (auto i=0;i<thread.size();++i){
        auto nid=thread[i].node;
        for (auto &l:trg.get_nodeview(nid).prev()){
            if (nodepos_in_reads.count(l.support().id)==0) nodepos_in_reads[l.support().id]={};
            else if (nodepos_in_reads[l.support().id].back().first==i){
                if ( nodepos_in_reads[l.support().id].back().second>l.support().index)
                    nodepos_in_reads[l.support().id].back().second=l.support().index;
                continue;
            }

            nodepos_in_reads[l.support().id].emplace_back(i,l.support().index);
        }
        for (auto &l:trg.get_nodeview(nid).next()){
            if (nodepos_in_reads.count(l.support().id)==0) nodepos_in_reads[l.support().id]={};
            else if (nodepos_in_reads[l.support().id].back().first==i){
                if ( nodepos_in_reads[l.support().id].back().second>l.support().index)
                    nodepos_in_reads[l.support().id].back().second=-l.support().index; //indexes are negative if the node is bw in thread
                continue;
            }

            nodepos_in_reads[l.support().id].emplace_back(i,l.support().index);
        }
    }
    //First check: chimeric junctions -> i.e links that are crossed by very few reads
    std::vector<uint64_t> max_read_coverage(thread.size()-1);//coverage across each transition of the thread
    for (auto &rnp:nodepos_in_reads){
        for (auto i=rnp.second.front().first;i<rnp.second.back().first;++i) ++max_read_coverage[i];
    }
    //for (auto &mrc:max_read_coverage) std::cout<<" "<<mrc;
    //std::cout<<std::endl;
    //TODO: we could "rescue" chimeras by splitting at the position of unsupported linkage
    for (auto i=5;i<max_read_coverage.size()-5;++i) if (max_read_coverage[i]<3) return {};

    //std::vector<uint64_t> happy_fw(thread.size());
    //std::vector<uint64_t> happy_bw(thread.size());
    //std::vector<uint64_t> unhappy_fw(thread.size());
    //std::vector<uint64_t> unhappy_bw(thread.size());
    //std::vector<uint64_t> disconnected_fw(thread.size());
    //std::vector<uint64_t> disconnected_bw(thread.size());
    std::vector<uint64_t> happy(thread.size());
    std::vector<uint64_t> unhappy(thread.size());
    std::vector<uint64_t> disconnected(thread.size());

    //TODO: this is a very fast approximation to the thing, but not as bad as it could be, could be better if considering ends of read
    for (auto &rnps:nodepos_in_reads){
        for (auto i=0;i<rnps.second.size();++i){
            //assess each node's happiness in this read
            Happiness fh=Happiness::unknown;
            Happiness bh=Happiness::unknown;

            if (i!=0){
                if ((rnps.second[i-1].second>0) == (rnps.second[i].second>0) and rnps.second[i-1].second<rnps.second[i].second)
                    bh=Happiness::happy;
                else
                    bh=Happiness::unhappy;
            }

            if (i!=rnps.second.size()-1){
                if ((rnps.second[i].second>0) == (rnps.second[i+1].second>0) and rnps.second[i].second<rnps.second[i+1].second)
                    fh=Happiness::happy;
                else
                    fh=Happiness::unhappy;
            }
            if (fh==Happiness::unhappy or bh==Happiness::unhappy) ++unhappy[rnps.second[i].first];
            else if (fh==Happiness::happy or bh==Happiness::happy) ++happy[rnps.second[i].first];
            else ++disconnected[rnps.second[i].first];
        }
    }
    //now just copy the components that are still happy
    std::vector<NodePosition> happy_nodes;
    happy_nodes.reserve(thread.size());
    for (auto i=0;i<thread.size();++i){
        if (unhappy[i]>max_unhappy) continue;
        if ((happy[i]+unhappy[i]+disconnected[i])*disconnection_rate<disconnected[i]) continue;
        happy_nodes.push_back(thread[i]);
    }
    return happy_nodes;
}

void make_all_threads_happy(LongReadsRecruiter & lrr, DistanceGraph &trg, int max_unhappy, float disconnection_rate){
    std::atomic<uint64_t> ptc(0);
#pragma omp parallel
    {
        uint64_t tptc=0;
#pragma omp for schedule(dynamic, 10000)
        for (auto i = 0; i < lrr.read_threads.size(); ++i) {
            if (++tptc==10000) {
                if (++ptc % 10 == 0) sdglib::OutputLog() << ptc*10000 << " threads processed" << std::endl;
                tptc=0;
            }
            if (lrr.read_threads[i].size() > 20) {
                lrr.read_threads[i] = make_thread_happy(lrr.read_threads[i], trg, max_unhappy, disconnection_rate);
            }
            else lrr.read_threads[i] = {};
        }
    }
}

std::unordered_map<sgNodeID_t , int64_t > sort_cc(const DistanceGraph& dg, std::unordered_set<sgNodeID_t> cc){
    // Next nodes to process
    std::vector<sgNodeID_t > next_nodes;
    // Map with the nodes positions
    std::unordered_map<sgNodeID_t , int64_t > node_pos;

    // check no node on the CC appears in both directions.
    for (const auto&x: cc){
        if (cc.find(-x) != cc.end()){
            return {};
        }
    }

    // now start from each node without prev as 0
    // Nodes with no prev are the inital nodes in the partial order
    for (const auto& x: cc){
        if (dg.get_nodeview(x).prev().empty()){
            next_nodes.push_back(x);
        }
        node_pos[x]=0;
    }


    // Propagate position FW
    auto pass=0;
    while (!next_nodes.empty()) {
        std::vector<sgNodeID_t > new_next_nodes;
        for (const auto& nid: next_nodes){
            for (const auto& l: dg.get_nodeview(nid).next()){
                // TODO: does this work like this --> if nid node ends after the linked node the next node is moved fw and the linked node is added to the list to place
                if (node_pos[nid]+dg.get_nodeview(nid).size()+l.distance()>node_pos[l.node().node_id()]){ // If the node position lands fw of the next node the next node position is pushed FW
                    node_pos[l.node().node_id()] = node_pos[nid]+dg.get_nodeview(nid).size()+l.distance();
                    new_next_nodes.push_back(l.node().node_id());
                }
            }
        }
        next_nodes=new_next_nodes;
        if (++pass>cc.size()){
            return {};
        };
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
    pass=0;
    while (updated){
        updated=false;
        //sort node_pos by the second element
        // TODO: This sorted_node_pos could be copies once to the vector outside the loop and then work with the vector until return
        std::vector<std::pair<sgNodeID_t , int64_t>> sorted_node_pos;
        std::copy(node_pos.begin(), node_pos.end(), std::back_inserter< std::vector<std::pair<sgNodeID_t , int64_t>>>(sorted_node_pos));
        std::sort(sorted_node_pos.begin(), sorted_node_pos.end(), [](std::pair<sgNodeID_t , int64_t> elem1 ,std::pair<sgNodeID_t , int64_t> elem2){return elem1.second < elem2.second;});

        for (const auto& mp: sorted_node_pos){
            auto nid = mp.first;
            auto p = node_pos[mp.first];
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
        if (++pass>cc.size()){
            return {};
        };
    }

    // move every node up to its limit bw
    updated=true;
    pass=0;
    while(updated){
        updated=false;
        // TODO: This sorted_node_pos could be copies once to the vector outside the loop and then work with the vector until return
        std::vector<std::pair<sgNodeID_t , int64_t>> sorted_node_pos;
        std::copy(node_pos.begin(), node_pos.end(), std::back_inserter< std::vector<std::pair<sgNodeID_t , int64_t>>>(sorted_node_pos));
        std::sort(sorted_node_pos.begin(), sorted_node_pos.end(), [](std::pair<sgNodeID_t , int64_t> elem1 ,std::pair<sgNodeID_t , int64_t> elem2){return elem1.second < elem2.second;});

        for (const auto& mp: sorted_node_pos){
            auto nid = mp.first;
            auto p = node_pos[mp.first];
            auto nid_node_size = dg.get_nodeview(nid).size();

            if (!dg.get_nodeview(nid).prev().empty()){
                // get min distance between the nid and the others, means the next node
                int64_t np = 0;
                for (const auto& l: dg.get_nodeview(nid).prev()){
                    auto linked_node_id = l.node().node_id();
                    if (node_pos[linked_node_id]+l.distance()+l.node().size() > np)
                        np = node_pos[linked_node_id]+l.distance()+l.node().size();
                }
                if (np!=node_pos[nid]){
                    node_pos[nid]=np;
                    updated=true;
                }
            }
        }
        if (++pass>cc.size()){
            return {};
        };
    }
    return node_pos;
}

bool pop_node(DistanceGraph& dg, sgNodeID_t node_id, uint64_t read){
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
        dg.remove_link(-node_id, ln[0].node().node_id(), ln[0].distance(), ln[0].support());
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

HappyPosition NodeAdjacencies::happy_to_add(float used_perc) {
    bool h_fw_prevs (prevs.size()*used_perc<=fw_prevs);
    bool h_bw_prevs (prevs.size()*used_perc<=bw_prevs);
    bool h_fw_nexts (nexts.size()*used_perc<=fw_nexts);
    bool h_bw_nexts (nexts.size()*used_perc<=bw_nexts);
    if (h_fw_nexts and h_bw_prevs) return Nowhere;
    if (h_bw_nexts and h_fw_prevs) return Nowhere;
    if (h_fw_nexts){
        if (fw_prevs) return MiddleFW;
        else return FrontFW;
    }
    if (h_bw_nexts){
        if (bw_prevs) return MiddleBW;
        else return FrontBW;
    }
    if (h_fw_prevs){
        if (fw_nexts) return MiddleFW;
        else return BackFW;
    }
    if (h_bw_prevs){
        if (bw_nexts) return MiddleBW;
        else return BackBW;
    }
    return Nowhere;
}

void NodeAdjacencies::mark_used(const sgNodeID_t &nid) {
    if (prevs.count(nid)) ++fw_prevs;
    if (prevs.count(-nid)) ++bw_prevs;
    if (nexts.count(nid)) ++fw_nexts;
    if (nexts.count(-nid)) ++bw_nexts;
}

void NodeAdjacencies::mark_unused(const sgNodeID_t &nid) {
    if (prevs.count(nid)) --fw_prevs;
    if (prevs.count(-nid)) --bw_prevs;
    if (nexts.count(nid)) --fw_nexts;
    if (nexts.count(-nid)) --bw_nexts;
}

std::pair<int64_t, int64_t> NodeAdjacencies::find_happy_place(const HappyInsertionSorter &sorter) {
    int64_t max_prev=INT64_MIN,min_next=INT64_MAX;
    //First check if the node is fw or backward
    if ((fw_prevs or fw_nexts) and (bw_nexts or bw_prevs)) return {max_prev,min_next}; //mixed orientations, unhappy
    //fordard
    if ((fw_prevs or fw_nexts)) {
        for (auto &nid:prevs) {
            auto p = sorter.get_node_position(nid);
            if (p == 0) continue;
            max_prev = std::max(max_prev, p);
        }
        for (auto &nid:nexts) {
            auto p = sorter.get_node_position(nid);
            if (p == 0) continue;
            min_next = std::min(min_next, p);
        }
    }
    else {
        for (auto &nid:nexts) {
            auto p = sorter.get_node_position(nid);
            if (p == 0) continue;
            if (llabs(p)>llabs(max_prev)) max_prev = p;
        }
        for (auto &nid:prevs) {
            auto p = sorter.get_node_position(nid);
            if (p == 0) continue;
            if (llabs(p)<llabs(min_next)) min_next = p;
        }
    }
    return {max_prev,min_next};
}

void HappyInsertionSorter::compute_adjacencies(int min_links,int radius) {
    //count all connections only once, along threads.
    std::unordered_map<std::pair<sgNodeID_t, sgNodeID_t>,uint64_t> conns;
    for (auto &nv:rtg.get_all_nodeviews(false,false))
        adjacencies[nv.node_id()]=NodeAdjacencies();
    for (auto &tinf:rtg.thread_info){
        auto thread=rtg.get_thread(tinf.first);
        for (auto i1=0;i1<thread.size();++i1){
            auto e1=-thread[i1].node;
            for (auto i2=i1+1;i2<thread.size() and i2<i1+radius;++i2) {
                auto e2=thread[i2].node;
                if (e1<e2) ++conns[{e1,e2}];
                else ++conns[{e2,e1}];
            }
        }

    }
    //add every connection that passes the min_links to both nodes's appropriate sets
    for (auto &it:conns){
        if (it.second>=min_links) {
            if (it.first.first<0) adjacencies[-it.first.first].nexts.insert(it.first.second);
            else adjacencies[it.first.first].prevs.insert(-it.first.second);
            if (it.first.second<0) adjacencies[-it.first.second].nexts.insert(it.first.first);
            else adjacencies[it.first.second].prevs.insert(-it.first.first);
        }
    }
}

NodeAdjacencies & HappyInsertionSorter::get_node_adjacencies(sgNodeID_t nid) {
    if (adjacencies.count(llabs(nid))==0) throw std::runtime_error("node is disconnected");
    return adjacencies[llabs(nid)];
}

void HappyInsertionSorter::add_node(sgNodeID_t nid) {
    if (node_positions.count(nid) or node_positions.count(-nid)) throw std::runtime_error("node "+std::to_string(nid)+" is already in the order!");
    node_positions[llabs(nid)]=0;//TODO: positions are not being used right now
    candidates.erase(nid);
    candidates.erase(-nid);
    if (nid>0){
        for (auto &p:adjacencies[nid].prevs) if (node_positions.count(p)==0 and node_positions.count(-p)==0) candidates.insert(p);
        for (auto &n:adjacencies[nid].nexts) if (node_positions.count(n)==0 and node_positions.count(-n)==0) candidates.insert(n);
    }
    else {
        for (auto &p:adjacencies[-nid].prevs) if (node_positions.count(p)==0 and node_positions.count(-p)==0) candidates.insert(-p);
        for (auto &n:adjacencies[-nid].nexts) if (node_positions.count(n)==0 and node_positions.count(-n)==0) candidates.insert(-n);
    }
    for (auto &p:adjacencies[llabs(nid)].prevs) adjacencies[llabs(p)].mark_used(nid);
    for (auto &n:adjacencies[llabs(nid)].nexts) adjacencies[llabs(n)].mark_used(nid);
}

void HappyInsertionSorter::remove_node_from_everywhere(sgNodeID_t nid) {
    candidates.erase(nid);
    candidates.erase(-nid);
    if (node_positions.count(nid)) {
        node_positions.erase(nid);
        for (auto &p:adjacencies[nid].prevs) adjacencies[llabs(p)].mark_unused(nid);
    }
    if (node_positions.count(-nid)) {
        node_positions.erase(-nid);
        for (auto &p:adjacencies[nid].prevs) adjacencies[llabs(p)].mark_unused(-nid);
    }
}

int64_t HappyInsertionSorter::get_node_position(sgNodeID_t nid) const {
    auto it=node_positions.find(nid);
    if (it!=node_positions.end()) return it->second;
    it=node_positions.find(-nid);
    if (it!=node_positions.end()) return -it->second;
    return 0;
}

bool HappyInsertionSorter::insert_node(sgNodeID_t nid, float used_perc, bool solve_floating_by_rtg) {
    if (node_positions.empty()) {
        if (adjacencies.count(llabs(nid))==0) return false;
        add_node(nid);
        if (nid>0) node_positions[nid]=1;
        else node_positions[-nid]=-1;
        return true;
    }
    nid=llabs(nid);
    if (node_positions.count(nid) or node_positions.count(-nid)) return false;
    if (adjacencies[nid].happy_to_add(used_perc)==HappyPosition::Nowhere) return false;
    auto place=adjacencies[nid].find_happy_place(*this);
    int64_t np=0;
    if (place.first==INT64_MIN and place.second==INT64_MAX) {
        return false;
    }
    else if (place.first==INT64_MIN and abs(place.second)==1) {
        np=place.second;
    }
    else if (llabs(place.first)==node_positions.size() and place.second==INT64_MAX) {
        if (place.first<0) np=place.first-1;
        else np=place.first+1;
    }
    else if (std::signbit(place.first)!=std::signbit(place.second) or llabs(place.first)>llabs(place.second)) {
        return false;
    }
    else if (llabs(place.first)==llabs(place.second)-1) {//perfect place found!
        np=place.second;
    }
    //TODO: solve floating
    if (np==0) return false;
    for (auto &it:node_positions) {
        if (it.second>=llabs(np)) ++it.second;
        if (-it.second>=llabs(np)) --it.second;
    }
    if (np>0) add_node(nid);
    else add_node(-nid);
    node_positions[nid]=np;
    return true;
}

void HappyInsertionSorter::reset_positions() {
    candidates.clear();
    node_positions.clear();
    for (auto &a:adjacencies){
        a.second.fw_prevs=0;
        a.second.fw_nexts=0;
        a.second.bw_prevs=0;
        a.second.bw_nexts=0;
    }
}

TheGreedySorter::TheGreedySorter(const DistanceGraph& _trg_nt, sgNodeID_t founding_node):trg_nt(_trg_nt), dg(_trg_nt.sdg){


    // fill all nodes vector
    if (founding_node){
        for (const auto &nid:trg_nt.get_connected_component(founding_node,false))
            all_nodes.push_back(nid);
    }
    else {
        for (const auto &n: trg_nt.get_all_nodeviews(false, false))
            all_nodes.push_back(n.node_id());
    }

    for (const auto& nid:all_nodes){
        auto nv=trg_nt.get_nodeview(nid);
        for (const auto& lv: nv.next())
            all_reads.insert(lv.support().id);
        for (const auto& lv: nv.prev())
            all_reads.insert(lv.support().id);

    }

    std::cout << "finding read ends" << std::endl;
#pragma omp parallel for
    for (auto i=0;i<all_nodes.size();++i){
        const auto &nid=all_nodes[i];
        auto nv=trg_nt.get_nodeview(nid);
        // fill fw and bw rids support vectors
        std::unordered_set<sgNodeID_t > frids;
        for (const auto& x: nv.next())
            frids.insert(x.support().id);
        std::unordered_set<sgNodeID_t > brids;
        for (const auto& x: nv.prev())
            brids.insert(x.support().id);

        // Checks if this rid is at the end of the available support
        for (const auto& rid: frids){
            if (brids.find(rid)==brids.end()){
                #pragma omp critical
                {
                    if (read_ends.find(rid)==read_ends.end()){
                        read_ends[rid]={};
                }
                    read_ends[rid].push_back(nv.node_id());
                }
            }
        }
        // Checks if this rid is at the end of the available support
        for (const auto& rid: brids){
            if (frids.find(rid) == frids.end()) {
                #pragma omp critical
                {
                    if (read_ends.find(rid)==read_ends.end()) {
                        read_ends[rid]={};
                    }
                    read_ends[rid].push_back(-nv.node_id());
                }
            }
        }
    }

    std::cout << "TheGreedySorter created with " << all_nodes.size() << " nodes and "<< all_reads.size() <<" reads"<<std::endl;
}

void TheGreedySorter::update_read_nodes_in_order() {
    read_nodes_in_order.clear();
    std::set<uint64_t> rids;
    for (auto nvo:dg.get_all_nodeviews(false,false)){
        auto nv=trg_nt.get_nodeview(nvo.node_id());
        rids.clear();
        for (auto const &l:nv.next())
            rids.insert(l.support().id);
        for (auto const &l:nv.prev())
            rids.insert(l.support().id);
        for (auto rid:rids) read_nodes_in_order[rid]+=1;
    }
}

std::vector<uint64_t> TheGreedySorter::node_belonging_scores(int64_t nid) {
    auto nv=trg_nt.get_nodeview(nid);
    std::set<uint64_t> rids;
    for (auto const &l:nv.next())
        rids.insert(l.support().id);
    for (auto const &l:nv.prev())
        rids.insert(l.support().id);
    std::vector<uint64_t> scores;
    for (auto rid:rids)  scores.push_back(read_nodes_in_order[rid]);
    return scores;
}

std::vector<uint64_t > TheGreedySorter::rids_from_node(NodeView nv){
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

uint64_t TheGreedySorter::shared_reads(NodeView nv1, NodeView nv2){
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

bool TheGreedySorter::pop_node(sgNodeID_t node_id, uint64_t read){
    order_is_valid=false;
    return ::pop_node(dg,node_id,read);
}

bool TheGreedySorter::remove_node(sgNodeID_t node_id) {
    auto nv=dg.get_nodeview(node_id);
    bool popped=false;
    for (auto l:nv.next()){
        if (pop_node(node_id,l.support().id)) popped=true;
    }

    for (auto l:nv.prev()){
        if (pop_node(node_id,l.support().id)) popped=true;
    }
    used_nodes.erase(node_id);
    used_nodes.erase(-node_id);
    return popped;
}

void TheGreedySorter::start_from_read(uint64_t rid, int min_confirmation){
    for (auto nv:dg.get_all_nodeviews(false,false)){
        dg.disconnect_node(nv.node_id());
    }
    used_reads.clear();
    used_nodes.clear();
    order_is_valid=false;
    // only pick the safest nodes alongside a read (i.e. they have a minimum number of confirming other reads shared)

    if(read_ends[rid].empty()){
        std::cout<< "No reads to start from!." << std::endl;
        return;
    }

    auto nv = trg_nt.get_nodeview(read_ends[rid][0]);
    int links_added=0;
    while (true){

        std::vector<LinkView> vlf;
        for (const auto& x: nv.next()){
            if (x.support().id == rid){
                vlf.push_back(x);
            }
        }
        if (vlf.empty()) break;
        auto lf = vlf[0];

//        std::cout << "Adding link" << -nv.node_id() << "-->" << lf.node().node_id() << "(" << lf.distance() << " " << lf.support().id << ")" << std::endl;
        dg.add_link(-nv.node_id(), lf.node().node_id(), lf.distance(), lf.support());
        nv = lf.node();
        links_added++;
    }
    std::cout << links_added<< " links added" << std::endl;
    // Pop all nodes with no minim support confirmation bw or fw
    bool popped=true;
    int nodes_popped=0;
    while (popped) {
        popped=false;
        int nvs_analysed = 0;
        for (const auto& nn: dg.get_all_nodeviews(true, false)){
            nvs_analysed++;
            if (!nn.next().empty()){
                // shared reads checks support of the current nv id in the multigra that stores all links from the reads (that why the nv roundabout, nv from different graphs)
                if (shared_reads(trg_nt.get_nodeview(nn.node_id()), trg_nt.get_nodeview(nn.next()[0].node().node_id())) <= min_confirmation ) {
                    pop_node(nn.node_id(), rid);
                    popped=true;
                    nodes_popped++;
                }
            }
        }
        if (nvs_analysed<10){
            std::cout<< "Aborting since there are less than 10 nodes left" << std::endl;
            break;
        }
    }
    std::cout << nodes_popped << " nodes popped" << std::endl;
    for (const auto& nn: dg.get_all_nodeviews(false, false)){
        used_nodes.insert(nn.node_id());
    }
    used_reads.push_back(rid);
}

std::unordered_map<sgNodeID_t , int64_t > TheGreedySorter::sort_graph(){
    if (!order_is_valid){
        if (used_nodes.empty()){
            std::cout << "No nodes, no order" << std::endl;
            // No node, returns empty order
            return {};
        }
        auto cc =  dg.get_connected_component(*used_nodes.begin());
        order=sort_cc(dg, cc);
        order_is_valid=true;
    }
    return order;
}

std::pair<int, int> TheGreedySorter::evaluate_read(uint64_t rid, bool print_pos){
    // Walk through a read in the original graph, count how many nodes are in the sorted graph, then check their order
    std::unordered_map<sgNodeID_t , int64_t > current_node_pos;
    if (print_pos){
        current_node_pos = sort_graph();
    }

    auto nv = trg_nt.get_nodeview(read_ends[rid][0]);
    int used = 0;
    int unused = 0;
    std::unordered_set<sgNodeID_t > seen={llabs(nv.node_id())};
    while (true) {
        auto nnv = dg.get_nodeview(nv.node_id());
        if (!nnv.prev().empty() or !nnv.next().empty()){
            used++;
            if(print_pos){
                if (current_node_pos.find(nv.node_id()) != current_node_pos.end())
                    std::cout << current_node_pos[nv.node_id()] << std::endl;
                if (current_node_pos.find(-nv.node_id()) != current_node_pos.end())
                    std::cout << -current_node_pos[-nv.node_id()] << std::endl;
            }
        } else {
            unused++;
            if (print_pos) std::cout << "unplaced" << std::endl;
        }

        std::vector<LinkView > lf;
        for (const auto& x: nv.next()){
            if (x.support().id == rid){
                lf.push_back(x);
            }
        }
        if (lf.empty()) break;
        nv=lf[0].node();
        if (seen.find(nv.node_id()) != seen.end()){
            return std::make_pair(0,0);
        } else {
            seen.insert(nv.node_id());
        }
    }
    return std::make_pair(used, unused);
}

std::vector<int64_t> TheGreedySorter::thread_nodes(int64_t rid) {
    std::vector<int64_t > tn;
    auto nv = trg_nt.get_nodeview(rid>0 ? read_ends[rid][0]:read_ends[-rid][1]);
    tn.emplace_back(nv.node_id());
    for (int c=1;tn.size()==c;++c){
        for (const auto& x: nv.next()){
            if (x.support().id==abs(rid)) {
                nv=x.node();
                tn.emplace_back(nv.node_id());
                break;
            }
        }
    }
    return tn;
}

std::vector<int64_t > TheGreedySorter::thread_node_positions(int64_t rid){
    std::vector<int64_t > np;
    auto nv = trg_nt.get_nodeview(rid>0 ? read_ends[rid][0]:read_ends[-rid][1]);
    auto current_node_pos = sort_graph();

    while (true) {
        // in the original it's none, here i'm using the null node
        if (current_node_pos.find(nv.node_id()) != current_node_pos.end()){
            np.emplace_back(current_node_pos[nv.node_id()]);
        }
        else if (current_node_pos.find(-nv.node_id()) != current_node_pos.end()) {
            np.emplace_back(-current_node_pos[-nv.node_id()]);
        }
        else{
            np.emplace_back(-10000000000);
        }

        std::vector<LinkView > lf;
        for (const auto& x: nv.next()){
            if (x.support().id == abs(rid)){
                lf.push_back(x);
            }
        }
        if (lf.empty()) break;
        nv = lf[0].node();
    }
    return np;
}

std::vector<int32_t > TheGreedySorter::evaluate_read_nodeorder(uint64_t rid, bool print_pos){
    auto nv = trg_nt.get_nodeview(read_ends[rid][0]);
    auto current_node_pos = sort_graph();

    int32_t ordered = 0;
    int32_t unordered = 0;
    int32_t unused = 0;

    uint32_t last_pos = 0;
    std::unordered_set<sgNodeID_t > seen={llabs(nv.node_id())};
    // TODO: check here i'm using 0 as a None, this might be the cause of the node difference or more!!!!!! <-------
    while (true) {
        // in the original it's none, here i'm using the null node
        int32_t p=0;
        if (current_node_pos.find(nv.node_id()) != current_node_pos.end()){
            p = current_node_pos[nv.node_id()];
        }
        else if (current_node_pos.find(-nv.node_id()) != current_node_pos.end()) {
            p = -current_node_pos[-nv.node_id()];
        }
        else{
            p = 0;
        }

        if (p==0){
            unused++;
        } else if (last_pos ==0 or last_pos<p) {
            last_pos = p;
            ordered++;
        } else {
            std::cout << "Node "<< nv.node_id() <<" is unordered in read "<< rid << std::endl;
            unordered++;
        }

        std::vector<LinkView > lf;
        for (const auto& x: nv.next()){
            if (x.support().id == rid){
                lf.push_back(x);
            }
        }
        if (lf.empty()) break;
        nv = lf[0].node();
        if (seen.find(nv.node_id())!=seen.end()){ // TODO: check the abs or not nv.node_id() if bug or not
            return std::vector<int32_t > {0, 100, 0};
        } else {
            seen.insert(nv.node_id()); // TODO: check the abs or not nv.node_id() if bug or not
        }
    }
    return {ordered, unordered, unused};
}

void TheGreedySorter::add_read(uint64_t rid, int min_confirmation){
    order_is_valid=false;
    auto nv = trg_nt.get_nodeview(read_ends[rid][0]);
    int links_added=0;
    std::vector<sgNodeID_t > added_nodes;
    while (true) {
        added_nodes.push_back(nv.node_id());
        added_nodes.push_back(-nv.node_id());

        std::vector<LinkView > vlf;
        for (const auto& x: nv.next()){
            if (x.support().id == rid){
                vlf.push_back(x);
            }
        }
        if (vlf.empty()) break;
        LinkView lf=vlf[0];

        dg.add_link(-nv.node_id(), lf.node().node_id(), lf.distance(), lf.support());
        nv=lf.node();
        links_added++;
    }
    std::cout << links_added << " links added" << std::endl;
    // Pop all nodes that dont have support confirmations bw or fw
    bool popped = true;
    int nodes_popped = 0;
    while (popped) {
        popped = false;
        int nvs_analysed = 0;
        for (const auto& nid: added_nodes){
            nv = dg.get_nodeview(nid);
            std::vector<LinkView > lf;
            for (const auto& x: nv.next()){
                if (x.support().id == rid){
                    lf.push_back(x);
                }
            }
            if (lf.empty()) continue;
            nvs_analysed++;
            if(shared_reads(trg_nt.get_nodeview(nv.node_id()), trg_nt.get_nodeview(lf[0].node().node_id()))<=min_confirmation){
                pop_node(nv.node_id(), rid);
                popped=true;
                nodes_popped++;
            }
        }
        if (nvs_analysed<2){
            std::cout << "aborting since <2 nodes are left" << std::endl;
            break;
        }
    }
    std::cout << nodes_popped << " nodes popped" << std::endl;
    for (const auto& nv: dg.get_all_nodeviews(false, false)){
        used_nodes.insert(nv.node_id());
    }
    used_reads.push_back(rid);
}

//void TheGreedySorter::extend_solution(int min_support, int min_shared, int min_new){
//    bool finished=false;
//    for (auto step=0; step<1000;++step){
//        std::cout << "Starting expansio step --> " << step << std::endl;
//        auto s = sort_graph()[1];
//        std::cout << "Current sorted size: " << s[-1].second - s.begining().second << std::endl;
//        int32_t added = 0;
//
////        [[x[0] for x in s[:30]],[x[0] for x in s[-30:]]]
//        std::vector<std::map<sgNodeID_t , int64_t >> s_f30;
//        int cont=0;
//        for (const auto& x: s){
//            if (con)
//        }
//
//        for (const auto& x: )
//
//    }
//}

void TheGreedySorter::write_connected_nodes_graph(std::string filename){
    std::vector<sgNodeID_t > sn;

    for (const auto& nv: dg.get_all_nodeviews(false, false)){
        if (!dg.links[llabs(nv.node_id())].empty()){
            sn.push_back(nv.node_id());
        }
    }
    dg.write_to_gfa1(filename, sn);
}

std::pair<sgNodeID_t,sgNodeID_t> TheGreedySorter::get_thread_ends(int64_t rid){
    if (read_ends.count(llabs(rid))==0) return {0,0};
    auto &v=read_ends[llabs(rid)];
    return {v[0],v[1]};
}
