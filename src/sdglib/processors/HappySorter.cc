//
// Created by Bernardo Clavijo (EI) on 12/02/2021.
//

#include "HappySorter.hpp"

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
            if (l.source==-nid) {
                ++total;
                if (threads.count(rtg.thread_fw_in_node(l.support.id,nid) ? l.support.id : -l.support.id)) ++happy;
            }
        }
    }
    else if (not prev and next) {
        for (auto const &l:rtg.links[llabs(nid)]) {
            if (l.source==nid) {
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
    for (auto &nc:node_count){
        if (nc.second>=min_threads and node_happiness(nc.first,true,true,min_threads)>min_happiness) candidates.insert(nc.first);
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
    auto hs=HappySorter(rtg);
    hs.order=LocalOrder(rtg.order_nodes(nodes));
    hs.recruit_all_happy_threads(.1);
    while (last_order.size()<hs.order.size()){
        last_order=hs.order;
        nodes=hs.order.as_signed_nodes();
        for (auto c:hs.find_internal_candidates()) nodes.push_back(c);
        hs.order=LocalOrder(rtg.order_nodes(nodes));
        hs.recruit_all_happy_threads();
    }
    threads.clear();
    fw_open_threads.clear();
    bw_open_threads.clear();
    order=LocalOrder(rtg.order_nodes(nodes));
    recruit_all_happy_threads();
    return;
}

bool HappySorter::grow_fw(int min_threads, bool verbose) {
    if (verbose) std::cout<<std::endl<<"New fw grow round starting, starting with an order of "<<order.size()<<" nodes"<<std::endl;

    //STEP 1 - find and order fw candidates
    std::vector<sgNodeID_t> candidates;
    for (auto &c:find_fw_candidates(min_node_happiness,min_threads)) candidates.emplace_back(c);
    if (verbose) std::cout<<"found "<<candidates.size()<<" candidates forward, including "<<order_end_size<<" in the order end"<<std::endl;
    auto sorted_candidates=rtg.order_nodes(candidates);
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

    sorted_candidates=rtg.order_nodes(candidates);
    if (verbose) std::cout<<"sorted and internal candidates size: "<<sorted_candidates.size()<<std::endl;
    if (sorted_candidates.size()<20) {
        if (verbose) std::cout<<"aborting internal recruiting, order failed!"<<std::endl;
        return true;
    }

    //STEP 5 - a very crude merge again

    //TODO: this is stupidly innefficient!
    if (verbose) std::cout<<"Current order has "<<order.size()<<" nodes"<<std::endl;
    new_nodes.resize(old_nodes_size);
    new_nodes.insert(new_nodes.end(),sorted_candidates.begin(),sorted_candidates.end());
    auto new_order=LocalOrder(new_nodes);
    if(new_order.node_positions.size() > 0){
        order=new_order;
    } else {
        return false;
    }
    if (verbose) std::cout<<"New order has "<<order.size()<<" nodes"<<std::endl;

    // STEP 6 - recruit more threads again
    otc=threads.size();
    recruit_all_happy_threads();
    if (verbose) std::cout<<threads.size()-otc<<" new happy threads recruited"<<std::endl;

    return true;
}