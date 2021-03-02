//
// Created by Bernardo Clavijo (EI) on 12/02/2021.
//

#include <chrono>
#include <random>
#include "HappySorter.hpp"

std::vector<sgNodeID_t> rtg_place_order(const ReadThreadsGraph & rtg,std::vector<sgNodeID_t> &nodes){
    auto placed=rtg.place_nodes({},nodes);
    std::vector<sgNodeID_t> order;
    order.reserve(placed.size());
    for (const auto &np:placed) order.emplace_back(np.first);
    return order;
}

void HappySorter::reverse() {
    order=order.reverse();
    std::unordered_set<sgNodeID_t> new_threads,new_fthreads,new_bthreads;
    for (auto tid:threads) new_threads.insert(-tid);
    std::swap(threads,new_threads);

    for (auto tid:fw_open_threads) new_bthreads.insert(-tid);
    for (auto tid:bw_open_threads) new_fthreads.insert(-tid);
    std::swap(fw_open_threads,new_fthreads);
    std::swap(bw_open_threads,new_bthreads);

    //TODO: translate coordinates: end-pos-node_size. (end is last pos + last size).


}

float HappySorter::thread_happiness(int64_t tid,int min_nodes) const {
    if (min_nodes==-1) min_nodes=min_thread_nodes;
    int64_t first_ti=-1,last_ti=-1,first_oi=-1,last_oi=-1,shared_nodes=0;
    auto tnp=rtg.get_thread(tid);
    for (auto i=0;i<tnp.size();++i) {
        auto p=order.get_node_position(tnp[i].node);
        if (p>0){
            if (first_ti==-1){
                first_ti=i;
                first_oi=p;
            }
            last_ti=i;
            last_oi=p;
            ++shared_nodes;
        }
    }
    if (last_ti==-1 or shared_nodes<min_nodes) {
        //std::cout<<"Thread happiness for thread "<<tid<<" = 0, only "<<shared_nodes<<"/"<<min_nodes<<" requered shared nodes"<<std::endl;
        return 0;
    }
    //std::cout<<"Thread happiness for thread"<<tid<<" = "<<((float) shared_nodes)/(last_ti+1-first_ti)<<std::endl<<" based on "<<shared_nodes<<" shared nodes"<<std::endl;
    return ((float) shared_nodes)/(last_ti+1-first_ti);
}

float HappySorter::node_happiness(sgNodeID_t nid, bool prev, bool next,int min_threads) const {
    if (min_threads==-1) min_threads=min_node_threads;
    uint64_t happy=0,total=0;
    if (nid==0 or rtg.sdg.links.size()-1<llabs(nid)) return 0;
    if (not prev and not next) return 0;
    else if (prev and not next) {
        for (auto const &l:rtg.links[llabs(nid)]) {
            if (l.source==nid) {
                ++total;
                if (threads.count(rtg.thread_fw_in_node(l.support.id,nid) ? l.support.id : -l.support.id)) ++happy;
            }
        }
    }
    else if (not prev and next) {
        for (auto const &l:rtg.links[llabs(nid)]) {
            if (l.source==-nid) {
                ++total;
                if (threads.count(rtg.thread_fw_in_node(l.support.id,nid) ? l.support.id : -l.support.id)) ++happy;
            }
        }
    }
    else {
        for (auto const &tid:rtg.node_threads(nid,true)){
            ++total;
            if (threads.count(tid)) ++happy;
        }

    }
    if (total<min_threads) return 0;
    return ((float) happy)/total;
}

void HappySorter::recruit_all_happy_threads(float min_happiness, int min_nodes) {
    if (min_nodes==-1) min_nodes=min_thread_nodes;
    if (min_happiness==-1) min_happiness=min_thread_happiness;
    for (const auto &np: order.node_positions){
        const auto &nid=np.second>0 ? np.first:-np.first;
        for (auto &tid:rtg.node_threads(nid,true)){
            if (threads.count(tid)==0 and threads.count(-tid)==0 and thread_happiness(tid,min_nodes)>=min_happiness) {
                threads.insert(tid);
                fw_open_threads.insert(tid);
                bw_open_threads.insert(tid);
            }
        }
    }
}


std::unordered_set<sgNodeID_t> HappySorter::find_fw_candidates(float min_happiness, int min_threads, int end_size) const {
    if (min_threads==-1) min_threads=min_node_threads;
    if (min_happiness==-1) min_happiness=min_node_happiness;
    if (end_size==-1) end_size=order_end_size;
    std::unordered_map<sgNodeID_t,int> node_links;
    std::unordered_set<sgNodeID_t> candidates;


    for (auto tid:fw_open_threads){
        auto tnps=rtg.get_thread(tid);
        int64_t last_inside_node=-1;
        for (auto i=0;i<tnps.size();++i) {
            if (order.get_node_position(tnps[i].node)>order.size()-end_size) {
                last_inside_node=i;
                break;
            }
        }
        for (auto i=last_inside_node;i<tnps.size();++i) ++node_links[tnps[i].node];
    }

    for (auto &nl:node_links) if (nl.second>=min_threads and node_happiness(nl.first,true,false,min_threads)>=min_happiness) candidates.insert(nl.first);

    //remove candidates that appear in both directions (should be rare?), also remove candidates located before end_size
    std::set<sgNodeID_t> to_delete;
    for (auto c:candidates) {
        if ( candidates.count(-c)) to_delete.insert(llabs(c));
        else {
            auto p=order.get_node_position(c);
            if (p<0 or (p>0 and p<=order.size()-end_size)) to_delete.insert(c);
        }
    }
    for (auto cd:to_delete) {
        candidates.erase(cd);
        candidates.erase(-cd);
    }
    return candidates;
}


std::unordered_set<sgNodeID_t> HappySorter::find_bw_candidates(float min_happiness, int min_threads, int end_size) const {
    if (min_threads==-1) min_threads=min_node_threads;
    if (min_happiness==-1) min_happiness=min_node_happiness;
    if (end_size==-1) end_size=order_end_size;
    std::unordered_map<sgNodeID_t,int> node_links;
    std::unordered_set<sgNodeID_t> candidates;


    for (auto tid:bw_open_threads){
        auto tnps=rtg.get_thread(tid);
        int64_t last_inside_node=-1;
        for (int i=tnps.size()-1;i>=0;--i) {
            auto o=order.get_node_position(tnps[i].node);
            if (o>0 and o<end_size) {
                last_inside_node=i;
                break;
            }
        }
        for (int i=last_inside_node;i>=0;--i) ++node_links[tnps[i].node];
    }

    for (auto &nl:node_links) if (nl.second>=min_threads and node_happiness(nl.first,false,true,min_threads)>=min_happiness) candidates.insert(nl.first);

    //remove candidates that appear in both directions (should be rare?), also remove candidates located before end_size
    std::set<sgNodeID_t> to_delete;
    for (auto c:candidates) {
        if ( candidates.count(-c)) to_delete.insert(llabs(c));
        else {
            auto p=order.get_node_position(c);
            if (p<0 or (p>0 and p>=end_size)) to_delete.insert(c);
        }
    }
    for (auto cd:to_delete) {
        candidates.erase(cd);
        candidates.erase(-cd);
    }
    return candidates;
}

std::unordered_set<sgNodeID_t> HappySorter::find_internal_candidates(float min_happiness, int min_threads,
                                                                     int32_t first, int32_t last) const {
    if (min_threads==-1) min_threads=min_node_threads;
    if (min_happiness==-1) min_happiness=min_node_happiness;
    if (last>order.size()+1) last=order.size();
    if (first<1) first=1;

//    //Step #1, find first and last nodes on each thread that are in the order
    std::map<int64_t, std::pair<int32_t,int32_t>> thread_limits;
    auto snodes=order.as_signed_nodes();
    if (snodes.size()!=order.size()) return {};
    for (auto i=first;i<=last;++i) {
        auto &nid=snodes[i-1];
        for (auto ntp:rtg.node_threadpositions(nid)) {
            if (threads.count(ntp.first)==0) continue;
            auto &tl=thread_limits[ntp.first];
            if (tl.first==0) {
                tl = {ntp.second, ntp.second};
            }
            else {
                if (tl.first>ntp.second) tl.first=ntp.second;
                if (tl.second<ntp.second) tl.second=ntp.second;
            }
        }
    }
//    //Step #2, count nodes between first and last
    std::unordered_map<sgNodeID_t,int64_t> node_count;
    std::unordered_set<sgNodeID_t> candidates;
    for (auto &tl:thread_limits){
        auto tnps=rtg.get_thread(tl.first);
        for (auto i=tl.second.first-1;i<tl.second.second;++i){
            if (order.node_positions.count(llabs(tnps[i].node))==0) //skip nodes already in order
                ++node_count[tnps[i].node];
        }
    }
    //add candidates meeting criteria, if reverse present already delete reverse and don't insert.
    for (auto &nc:node_count){
        if (nc.second>=min_threads and node_happiness(nc.first,true,true,min_threads)>min_happiness) {
            if (candidates.count(-nc.first)) candidates.erase(-nc.first);
            else candidates.insert(nc.first);
        }

    }
    return candidates;
}

void HappySorter::close_internal_threads(int order_end, int thread_end) {
    std::set<uint64_t> to_close;
    for (auto tid:fw_open_threads){
        auto tnp=rtg.get_thread(tid);
        for (int i=tnp.size()-1;i>=0;--i){
            auto nid=tnp[i].node;
            auto p=order.get_node_position(nid);
            if (p>0) {
                if (tnp.size()-1-i<=thread_end or order.size()-p>order_end) to_close.insert(tid);
                break;
            }
        }
    }
    for (auto tid:to_close) fw_open_threads.erase(tid);
    to_close.clear();
    for (auto tid:bw_open_threads){
        auto tnp=rtg.get_thread(tid);
        for (int i=0;i<tnp.size();++i){
            auto nid=tnp[i].node;
            auto p=order.get_node_position(nid);
            if (p>0) {
                if (i<=thread_end or p-1>order_end) to_close.insert(tid);
                break;
            }
        }
    }
    for (auto tid:to_close) bw_open_threads.erase(tid);
}

void HappySorter::start_from_node(sgNodeID_t nid, int min_links) {
    std::unordered_map<sgNodeID_t,int> node_links;
    for (auto tid:rtg.node_threads(nid,true))
        for (auto &ntp:rtg.get_thread(tid))
            ++node_links[ntp.node];
    std::vector<sgNodeID_t> nodes;
    LocalOrder last_order;
    for (auto &nl:node_links) if (nl.second>=min_links) nodes.push_back(nl.first);
    auto hs=HappySorter(rtg,min_thread_happiness,min_node_happiness,min_thread_nodes,min_node_threads);
    //hs.order=LocalOrder(rtg.order_nodes(nodes));
    hs.order=LocalOrder(rtg_place_order(rtg,nodes));
    hs.recruit_all_happy_threads(.1);
    while (last_order.size()<hs.order.size()){
        last_order=hs.order;
        nodes=hs.order.as_signed_nodes();
        for (auto c:hs.find_internal_candidates()) nodes.push_back(c);
        //hs.order=LocalOrder(rtg.order_nodes(nodes));
        hs.order=LocalOrder(rtg_place_order(rtg,nodes));
        hs.recruit_all_happy_threads();
    }
    threads.clear();
    fw_open_threads.clear();
    bw_open_threads.clear();
    //order=LocalOrder(rtg.order_nodes(nodes));
    //order=LocalOrder(rtg_place_order(rtg,nodes));

    //auto placed=rtg.place_nodes({},nodes);
    order.node_positions[nid]=1;
    node_coordinates[nid]=0;
    std::unordered_set<sgNodeID_t> nodeset(nodes.begin(),nodes.end());
    add_placed_nodes(place_nodes(nodeset,false));
    recruit_all_happy_threads();
    close_internal_threads();
    return;
}

bool HappySorter::grow_fw(int min_threads, bool verbose) {
    if (verbose) std::cout<<std::endl<<"New fw grow round starting, starting with an order of "<<order.size()<<" nodes"<<std::endl;

    //STEP 1 - find and order fw candidates
    std::vector<sgNodeID_t> candidates;
    for (auto &c:find_fw_candidates(min_node_happiness,min_threads)) candidates.emplace_back(c);
    if (verbose) std::cout<<"found "<<candidates.size()<<" candidates forward, including "<<order_end_size<<" in the order end"<<std::endl;
    //auto sorted_candidates=rtg.order_nodes(candidates);
    auto sorted_candidates=rtg_place_order(rtg,candidates);
    if (verbose) std::cout<<"sorted candidates size: "<<sorted_candidates.size()<<std::endl;
    if (sorted_candidates.size()<20) {
        if (verbose) std::cout<<"aborting grow, not enough sorted candidates"<<std::endl;
        return false;
    }

    //STEP 2 - a very crude merge

    //TODO: this is stupidly innefficient!
    auto new_nodes=order.as_signed_nodes();
    new_nodes.resize(new_nodes.size()-order_end_size);
    auto old_nodes_size=new_nodes.size();
    new_nodes.insert(new_nodes.end(),sorted_candidates.begin(),sorted_candidates.end());
    order=LocalOrder(new_nodes);
    if (order.size()!=new_nodes.size()){
        std::cout<<"LocalOrder from old order and candidates has less nodes that its input!"<<std::endl;
    }
    if (verbose) std::cout<<"New order has "<<order.size()<<" nodes"<<std::endl;

    // STEP 3 - recruit more threads
    auto otc=threads.size();
    recruit_all_happy_threads();
    if (verbose) std::cout<<threads.size()-otc<<" new happy threads recruited"<<std::endl;

    // STEP 4 -
    auto internal_candidates=find_internal_candidates(.1,min_node_threads,old_nodes_size);
    if (verbose) std::cout<<"There are "<<internal_candidates.size()<<" internal candidates in the newly ordered region"<<std::endl;
    for (auto ic:internal_candidates) candidates.emplace_back(ic);

    //sorted_candidates=rtg.order_nodes(candidates);
    sorted_candidates=rtg_place_order(rtg,candidates);
    if (verbose) std::cout<<"sorted and internal candidates size: "<<sorted_candidates.size()<<std::endl;
    if (sorted_candidates.size()<20) {
        if (verbose) std::cout<<"aborting internal recruiting, order failed!"<<std::endl;
        return true;
    }

    //STEP 5 - a very crude merge again

    //TODO: this is stupidly innefficient!
    new_nodes.resize(old_nodes_size);
    new_nodes.insert(new_nodes.end(),sorted_candidates.begin(),sorted_candidates.end());
    order=LocalOrder(new_nodes);
    if (verbose) std::cout<<"New order has "<<order.size()<<" nodes"<<std::endl;

    // STEP 6 - recruit more threads again
    otc=threads.size();
    recruit_all_happy_threads();
    if (verbose) std::cout<<threads.size()-otc<<" new happy threads recruited"<<std::endl;

    return true;
}

std::map<int64_t, std::vector<std::pair<int64_t, sgNodeID_t>>> HappySorter::make_thread_nodepositions(const std::unordered_set<sgNodeID_t> & nodes) const{
    std::map<int64_t, std::vector<std::pair<int64_t, sgNodeID_t>>> thread_node_positions;

    //make the list of threads
    std::set<int64_t> tids;
    for (auto nid:nodes){
        auto ntids=rtg.node_threads(nid,true);
        for (auto tid:ntids) tids.insert(tid);
    }

    //TODO: count all fw/bw in node for each thread and discard those where nodes appear in more than one direciton.
    for (auto &tid:tids){
        auto &tnp=thread_node_positions[tid];
        auto s_nv=rtg.thread_start_nodeview(tid);
        while (node_coordinates.count(s_nv.node_id())==0 and nodes.count(s_nv.node_id())==0) {
            s_nv = rtg.next_in_thread(s_nv.node_id(),llabs(tid)).node();
        }
        int last_end_p=rtg.sdg.get_node_size(s_nv.node_id());
        LinkView ln(s_nv,0,Support());//invalid, only to allow the try block
        try {
            ln=rtg.next_in_thread(s_nv.node_id(),llabs(tid));
        }
        catch (const std::exception&) {
            continue;
        }
        tnp.emplace_back(0,s_nv.node_id());
        while (true){
            if (node_coordinates.count(ln.node().node_id()) or nodes.count(ln.node().node_id())) {
                tnp.emplace_back(last_end_p+ln.distance(),ln.node().node_id());
            }

            last_end_p+=ln.distance()+ln.node().size();

            try {
                ln=rtg.next_in_thread(ln.node().node_id(),llabs(tid));
            }
            catch (const std::exception&) {
                break;
            }
        }
    }
    for (auto it=thread_node_positions.begin();it!=thread_node_positions.end();) {
        if (it->second.empty()) it=thread_node_positions.erase(it);
        else ++it;
    }
    return thread_node_positions;

}


int64_t hs_place_node(const std::unordered_map<sgNodeID_t, int64_t> &node_positions, const std::map<sgNodeID_t, std::vector<std::pair<sgNodeID_t,int64_t>>> & node_distances, sgNodeID_t nid){
    if (node_distances.count(nid)==0) return INT64_MIN; //unplaced
    std::vector<int64_t> positions;
    for (const auto & nd:node_distances.at(nid)) {
        if (node_positions.count(nd.first))
            positions.emplace_back(node_positions.at(nd.first)+nd.second);
    }
    if (positions.empty()) return INT64_MIN;
    std::sort(positions.begin(),positions.end());
    return positions[positions.size()/2]; //Poor man's median hits back
};

void hs_update_npcomplete(std::map<sgNodeID_t, std::pair<bool,bool>> &np_complete,const std::unordered_map<sgNodeID_t, int64_t> &node_positions, const std::map<sgNodeID_t, std::vector<std::pair<sgNodeID_t,int64_t>>> & node_distances, const std::unordered_set<sgNodeID_t> &to_place){
    for (auto nid:to_place){
        if (node_distances.count(nid)==0) {
            np_complete[nid]={false,false};
            continue;
        }
        bool prevs=true,nexts=true;
        for (auto nd:node_distances.at(nid)) { //TODO: check prevs/nexts only if still true? saves lookups
            if (node_positions.count(nd.first)==0){
                if (nd.second>0) nexts= false;
                else prevs=false;
            }
        }
        np_complete[nid]={prevs,nexts};
    }
}

sgNodeID_t hs_most_connected_node(const std::unordered_map<sgNodeID_t, int64_t> &node_positions, const std::map<sgNodeID_t, std::vector<std::pair<sgNodeID_t,int64_t>>> & node_distances, const std::unordered_set<sgNodeID_t> &to_place){
    int most_connected=0;
    sgNodeID_t best_nid=0;
    for (auto nid:to_place){
        if (node_positions.count(nid)) continue;
        int prevs=0,nexts=0;
        for (auto nd:node_distances.at(nid)) {
            if (node_positions.count(nd.first)){
                if (nd.second>0) ++nexts;
                else ++prevs;
            }
        }
        if (nexts+prevs>most_connected){
            best_nid=nid;
            most_connected=nexts+prevs;
        }
    }
    return best_nid;
}

std::map<sgNodeID_t, std::vector<std::pair<sgNodeID_t,int64_t>>> hs_tnp_to_distances (const std::map<int64_t, std::vector<std::pair<int64_t, sgNodeID_t>>> &thread_nodepositions,const std::unordered_set<sgNodeID_t> &nodeset){
    std::map<sgNodeID_t, std::vector<std::pair<sgNodeID_t,int64_t>>> distances;

    for (const auto &tnp:thread_nodepositions) {
        sgNodeID_t last_node=0;
        int64_t last_pos=0;
        for (auto p:tnp.second){
            if (last_node!=0) {
                if (nodeset.count(p.second)) distances[p.second].emplace_back(last_node,p.first-last_pos);
                if (nodeset.count(last_node)) distances[last_node].emplace_back(p.second,last_pos-p.first);
            }
            last_node=p.second;
            last_pos=p.first;
        }
    }
    return distances;
}

//TODO: should we weight distances?
std::vector<std::pair<sgNodeID_t, int64_t>> HappySorter::place_nodes( const std::unordered_set<sgNodeID_t> &nodeset, bool verbose) const {

    if (verbose) std::cout<<"Place nodes starting"<<std::endl;

    //STEP 1: populate node distances for each node
    auto full_nodeset=nodeset;
    for (auto &pn:node_coordinates) full_nodeset.insert(pn.first);
    if (verbose) std::cout<<"Nodeset has "<<nodeset.size()<<" nodes, full nodeset has "<<full_nodeset.size()<<std::endl;
    auto tnp=make_thread_nodepositions(nodeset);
    if (verbose) std::cout<<"Created thread nodepositions from "<<tnp.size()<<" threads"<<std::endl;
    auto node_distances=hs_tnp_to_distances(tnp,nodeset);
    if (verbose) std::cout<<"Created distances for "<<node_distances.size()<<" nodes"<<std::endl;

    if (verbose) {
        for (auto &nd:node_distances) {
            std::cout<<"Distances for node "<<nd.first<<":";
            for (auto d:nd.second) std::cout<<" ( "<<d.first<<": "<<d.second<<")";
            std::cout<<std::endl;
        }
    }

    //STEP 2: place each node in the median of its placed distances.
    std::unordered_map<sgNodeID_t, int64_t> node_positions=node_coordinates;




    auto to_place=nodeset;
    std::map<sgNodeID_t, std::pair<bool,bool>> np_complete;
    if (verbose) std::cout<<"Entering placing loop"<<std::endl;
    while (not to_place.empty()){
        hs_update_npcomplete(np_complete,node_positions,node_distances,to_place);
        std::set<sgNodeID_t> placed;
        //first place with both prev and nexts full; if none, with prevs full; if none, with nexts full
        for (std::pair<bool,bool> cond:{std::make_pair(true,true),{true,false},{false,true}}) {

            for (auto nid:to_place) {
                if (np_complete[nid] == cond) {
                    auto p = hs_place_node(node_positions, node_distances, nid);
                    if (p != INT64_MIN) {
                        node_positions[nid] = p;
                        placed.insert(nid);
                    }
                }
            }
            if (not placed.empty()) break;
        }


        //if none placed, abort (for now), we can change np_complete to count unsatisfied and place the smallest
        if (placed.empty()) {
            //find node with most connections, place that one
            //std::cout<<"Aborting after placing "<<node_positions.size()<<" nodes, with "<<to_place.size()<<" nodes still to place, but no candidates"<<std::endl;
            auto nid=hs_most_connected_node(node_positions,node_distances,to_place);
            auto p = hs_place_node(node_positions, node_distances, nid);
            if (p != INT64_MIN) {
                node_positions[nid] = p;
                placed.insert(nid);
            }
            else {
                if (verbose) std::cout<<"Aborting after placing "<<node_positions.size()<<" nodes, with "<<to_place.size()<<" nodes still to place, but no candidates"<<std::endl;
                break;
            }
        }
        for (auto nid:placed) to_place.erase(nid);
    }

    std::vector<std::pair<sgNodeID_t,int64_t>> placements;
    placements.reserve(nodeset.size());
    for (auto &np:node_positions) if (nodeset.count(np.first)) placements.emplace_back(np.first,np.second);
    std::sort(placements.begin(),placements.end(),[](auto &a,auto &b){return a.second<b.second;});

    return placements;
}

bool HappySorter::add_placed_nodes(const std::vector<std::pair<sgNodeID_t, int64_t>> &placed_nodes,
                                   bool update_current) {
    auto old_node_count=order.node_positions.size();
    for (auto pn:placed_nodes) {
        if (update_current or node_coordinates.count(pn.first)==0) node_coordinates[pn.first]=pn.second;
    }
    //std::vector<std::pair<sgNodeID_t, int64_t>> all_nodes(node_coordinates.begin(),node_coordinates.end());
    std::vector<std::pair<sgNodeID_t, int64_t>> all_nodes;
    for (auto &nc:node_coordinates) all_nodes.emplace_back(nc.first,nc.second);
    std::sort(all_nodes.begin(),all_nodes.end(),[](auto &a,auto &b){return a.second<b.second;});
    order.node_positions.clear();
    for (auto i=0;i<all_nodes.size();++i) order.node_positions[llabs(all_nodes[i].first)]=all_nodes[i].first>0? i+1:-i-1;
    return order.node_positions.size()>old_node_count;
}

bool HappySorter::grow(int min_threads, float min_happiness, bool fw, bool bw, bool internal) {
    if (min_threads==-1) min_threads=min_node_threads;
    if (min_happiness==-1) min_happiness=min_node_happiness;
    auto old_size=node_coordinates.size();
    std::unordered_set<sgNodeID_t> candidates;
    if (internal) {
        add_placed_nodes(place_nodes(find_internal_candidates(min_happiness, min_threads), false));
        recruit_all_happy_threads();
        close_internal_threads();
    }
    if (fw) {
        add_placed_nodes(place_nodes(find_fw_candidates(min_happiness, min_threads), false));
        recruit_all_happy_threads();
        close_internal_threads();
    }
    if (bw) {
        add_placed_nodes(place_nodes(find_bw_candidates(min_happiness, min_threads), false));
        recruit_all_happy_threads();
        close_internal_threads();
    }
    return node_coordinates.size()>old_size;
}

//TODO: add verbose printing of first/last coordinate and option to grow only fw or bw.
//TODO2: add option to restrict internal candidates to the last added nodes+ some length at the end, and only run it full every few cycles (or never)
bool HappySorter::grow_loop(int min_threads, float min_happiness, int64_t steps, bool verbose) {
    bool grown=false;
    int64_t step=0;
    while (++step<steps) {
        if (grow(min_threads, min_happiness)) grown=true;
        else break;
        if (verbose and step%10==0) {
            std::pair<sgNodeID_t,int64_t> cmax={0,0},cmin={0,0};
            for (const auto &nc:node_coordinates) {
                if (nc.second<cmin.second) cmin=nc;
                if (nc.second>cmax.second) cmax=nc;
            }
            std::cout<<"After "<<step<<" grow rounds: "<<node_coordinates.size()<<" anchors span "<<cmax.second-cmin.second<<" bp, "<<cmin.first<<" @ "<<cmin.second<<" : "<<cmax.first<<" @ "<<cmax.second<<std::endl;
        }
    }
    return grown;
}

void HappySorterRunner::run(int64_t min_starting_nodes, float max_starting_used, int64_t min_final_nodes, int64_t max_steps, int64_t max_orders) {
    sdglib::OutputLog()<<"Starting HappySorterRunner run..."<<std::endl;
    auto nvs=rtg.get_all_nodeviews(false,false);
    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();

    shuffle (nvs.begin(), nvs.end(), std::default_random_engine(seed));

#pragma omp parallel for schedule(dynamic,100)
    for (auto snvi=nvs.cbegin();snvi<nvs.cend();++snvi){
        auto &snv=*snvi;
        if (orders.size()>=max_orders) continue;
        try {
#pragma omp critical(node_sorted)
            auto ns = node_sorted[snv.node_id()];
            if (node_sorted[snv.node_id()]) continue;
            auto hs = HappySorter(rtg, min_thread_happiness, min_node_happiness, min_thread_nodes, min_node_threads,
                                  order_end_size);
            hs.start_from_node(snv.node_id());
            if (hs.order.size() < min_starting_nodes) continue;
            int64_t used = 0;
#pragma omp critical(node_sorted)
            for (auto &nid:hs.order.as_signed_nodes()) if (node_sorted[llabs(nid)]) ++used;
            if (hs.order.size() * max_starting_used < used) continue;
            std::cout << "Thread #" << omp_get_thread_num() << " growing order from node " << snv.node_id() << " ..."
                      << std::endl;
            hs.grow_loop(-1, -1, max_steps);
            if (hs.order.size() < min_final_nodes) continue;
            //sorters[snv.node_id()]=hs;
#pragma omp critical(node_sorted)
            {
                if (orders.size()<max_orders and not node_sorted[snv.node_id()]) { //check again to mantain initial condition
                    std::cout << "Thread #" << omp_get_thread_num() << " adding order from node " << snv.node_id()
                              << " with " << hs.order.size() << " nodes" << std::endl;
                    orders[snv.node_id()] = hs.order;
                    for (auto &nid:hs.order.as_signed_nodes()) {
                        node_sorted[llabs(nid)] = true;
                        if (nid > 0) node_orders[nid].emplace_back(snv.node_id());
                        else node_orders[-nid].emplace_back(-snv.node_id());
                    }
                    if (orders.size() % 10 == 0) {
                        sdglib::OutputLog() << orders.size() << " orders created, "
                                            << std::count(node_sorted.begin(), node_sorted.end(), true)
                                            << " nodes sorted" << std::endl;
                    }
                }
            }
        }
        catch (std::exception &e) {
            std::cout<<"Exception caught when trying to process node "<<snv.node_id()<<": "<< typeid(e).name() <<": "<<e.what()<<", skipping"<<std::endl;
        }

    }
    sdglib::OutputLog()<<"HappySorterRunner run finished!!!"<<std::endl;
}

void HappySorterRunner::dump(std::string filename) {
    std::ofstream of(filename);
    of.write((char *)&min_thread_happiness,sizeof(min_thread_happiness));
    of.write((char *)&min_thread_nodes,sizeof(min_thread_nodes));
    of.write((char *)&min_node_happiness,sizeof(min_node_happiness));
    of.write((char *)&min_node_threads,sizeof(min_node_threads));
    of.write((char *)&order_end_size,sizeof(order_end_size));
    uint64_t sval=orders.size();
    of.write((char *) &sval,sizeof(sval));
    for (auto &o:orders){
        of.write((char *) &o.first,sizeof(o.first));
        auto nodes=o.second.as_signed_nodes();
        sdglib::write_flat_vector(of,nodes);
    }
    sval=node_orders.size();
    of.write((char *) &sval,sizeof(sval));
    for (auto &o:node_orders){
        of.write((char *) &o.first,sizeof(o.first));
        sdglib::write_flat_vector(of,o.second);
    }
    //bool vector can't be written with the usual functions due to packing, expand it to uint8_t
    std::vector<uint8_t> nsv;
    nsv.reserve(node_sorted.size());
    for (const auto &ns:node_sorted) nsv.push_back(ns);
    sdglib::write_flat_vector(of,nsv);

}

void HappySorterRunner::load(std::string filename) {
    std::ifstream ifs(filename);
    ifs.read((char *)&min_thread_happiness,sizeof(min_thread_happiness));
    ifs.read((char *)&min_thread_nodes,sizeof(min_thread_nodes));
    ifs.read((char *)&min_node_happiness,sizeof(min_node_happiness));
    ifs.read((char *)&min_node_threads,sizeof(min_node_threads));
    ifs.read((char *)&order_end_size,sizeof(order_end_size));
    uint64_t sval;
    ifs.read((char *) &sval,sizeof(sval));
    sgNodeID_t nid;
    std::vector<sgNodeID_t> nodes;
    orders.clear();
    while (sval--){
        ifs.read((char *) &nid,sizeof(nid));
        sdglib::read_flat_vector(ifs,nodes);
        orders.emplace(nid,nodes);
    }

    ifs.read((char *) &sval,sizeof(sval));
    node_orders.clear();
    while (sval--){
        ifs.read((char *) &nid,sizeof(nid));
        sdglib::read_flat_vector(ifs,nodes);
        node_orders.emplace(nid,nodes);
    }

    //bool vector can't be written with the usual functions due to packing, read expanded uint8_t, then insert
    std::vector<uint8_t> nsv;
    sdglib::read_flat_vector(ifs,nsv);
    node_sorted.reserve(nsv.size());
    for (const auto &ns:nsv) node_sorted.push_back(ns);


}