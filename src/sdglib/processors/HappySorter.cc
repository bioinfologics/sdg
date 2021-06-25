//
// Created by Bernardo Clavijo (EI) on 12/02/2021.
//

#include <chrono>
#include <random>
#include "HappySorter.hpp"
#include <deque>



std::vector<sgNodeID_t> LocalOrder::as_signed_nodes() const {
    std::vector<sgNodeID_t> nodes(node_positions.size());
    for (auto &np:node_positions){
        if (llabs(np.second)>nodes.size()) return {};
        nodes[llabs(np.second)-1]=(np.second>0 ? np.first:-np.first);
    }
    return nodes;
}

std::vector<NodePosition> LocalOrder::as_thread(const DistanceGraph &dg) const {
    std::vector<NodePosition> thread;
    if (node_coordinates.empty()) return {};
    thread.reserve(node_coordinates.size());
    //first place only starts and sort
    for (auto &nc:node_coordinates)  thread.emplace_back(nc.first,nc.second,0);
    std::sort(thread.begin(),thread.end(),[](auto a, auto b){return a.start<b.start;});
    //now update to be 0-based and fill ends
    auto offset=thread.front().start;
    for (auto &np:thread) {
        np.start-=offset;
        np.end=np.start+dg.sdg.get_node_size(np.node);
    }
    return thread;
}

int64_t LocalOrder::get_node_position(sgNodeID_t nid) const {
    auto it=node_positions.find(llabs(nid));
    if (it==node_positions.end()) return 0;
    return nid>0 ? it->second:-it->second;
}
//TODO:implement the choosing a node its neighbours not in enough orders yet, up to a certain coverage, single threaded at first

LocalOrder::LocalOrder(const std::vector<sgNodeID_t> &nodes) {
    for (int64_t i=0;i<nodes.size();++i) {
        auto &p = node_positions[llabs(nodes[i])];
        if (p!=0) {
            std::cout<<"ERROR, avoiding creation of LocalOrder with repeated node "<<llabs(nodes[i])<<" at "<<p<<" and "<< (nodes[i] > 0 ? i + 1 : -i - 1)<<std::endl;
            node_positions.clear();
            return;
        }
        p = nodes[i] > 0 ? i + 1 : -i - 1;
    }
}

LocalOrder LocalOrder::merge(const LocalOrder &other,int max_overhang,float min_shared_perc,int min_shared,float max_disordered_perc) const{
    //first figure out direction of merge
    int64_t first_seen_at=-1, first_seen_other, last_seen_at, last_seen_other;
    bool last_disordered=false;//keep track of last disordered, we won't merge if first or last are disordered
    uint64_t ordered=0,disordered=0;
    auto nids=as_signed_nodes();
    bool same_direction;
    std::set<sgNodeID_t> disordered_nodes;
    for (int64_t i=0;i<nids.size();++i){
        auto p=other.node_positions.find(llabs(nids[i]));
        if (p==other.node_positions.end()) continue;
        auto op=nids[i]>0 ? p->second: -p->second;
        if (first_seen_at==-1){
            first_seen_at=i;
            first_seen_other=op;
            last_seen_at=i;
            last_seen_other=op;
            ++ordered;
            same_direction=( op>0 );
        }
        else {
            last_disordered=false;
            if (same_direction) {//same sign, both orders are going in the same direction
                if ( op>0 and llabs(op) > llabs(last_seen_other)) {
                    ++ordered;
                    last_seen_at=i;
                    last_seen_other=op;
                }
                else {
                    ++disordered;
                    last_disordered=true;
                }
            }
            else {
                if ( op<0 and llabs(op) < llabs(last_seen_other)) {
                    ++ordered;
                    last_seen_at=i;
                    last_seen_other=op;
                }
                else {
                    ++disordered;
                    disordered_nodes.insert(llabs(nids[i]));
                    last_disordered=true;
                }
            }
        }
    }

    if (first_seen_at==-1) return LocalOrder();
    ++first_seen_at;
    ++last_seen_at;

    int64_t this_left=first_seen_at-1,other_left=(same_direction ? first_seen_other-1:other.node_positions.size()+first_seen_other);
    int64_t this_right=node_positions.size()-last_seen_at,other_right=(same_direction ? other.node_positions.size()-last_seen_other : -last_seen_other-1 );

    auto largest_left=std::max(this_left,other_left);
    auto largest_right=std::max(this_right,other_right);
    float shared_perc=(float)ordered/std::max(last_seen_at+1-first_seen_at,(int64_t) (llabs(last_seen_other)+1-llabs(first_seen_other)));
    float disordered_perc=(float)disordered/llabs(std::min(last_seen_at+1-first_seen_at,(int64_t) (llabs(last_seen_other)+1-llabs(first_seen_other))));

//    sdglib::OutputLog()<<"First seen at: "<<first_seen_at<<std::endl;
//    sdglib::OutputLog()<<"First seen other: "<<first_seen_other<<std::endl;
//    sdglib::OutputLog()<<"Last seen at: "<<last_seen_at<<std::endl;
//    sdglib::OutputLog()<<"Last seen other: "<<last_seen_other<<std::endl;
//    sdglib::OutputLog()<<"Ordered: "<<ordered<<" disordered: "<<disordered<<std::endl;
//    sdglib::OutputLog()<<"This left:"<<this_left<<", other left:"<<other_left<<std::endl;
//    sdglib::OutputLog()<<"This right:"<<this_right<<", other right:"<<other_right<<std::endl;
//    sdglib::OutputLog()<<"Shared perc: "<<shared_perc<<std::endl;
//    sdglib::OutputLog()<<"Disordered perc: "<<disordered_perc<<std::endl;

    if (ordered<min_shared or shared_perc<min_shared_perc or disordered_perc>max_disordered_perc or last_disordered or std::min(this_right,other_right)>max_overhang or std::min(this_left,other_left)>max_overhang) return LocalOrder();

    /**
     * At this point: this_left/other left is tip sizes BEFORE this order, right is AFTER.
     * All coordinates are in 1 : size, with negatives representing revrse complement vs. the corresponding node IN THIS ORDER.
     * First/last seen are INCLUDED boundaries of the overlapping segment.
     */


    //reverse other order if needed (easier to just write one fw-merging algorithm)
    auto other_fw=other;
    if (not same_direction) {
        other_fw=other.reverse();
        first_seen_other=other.node_positions.size()+1+first_seen_other;
        last_seen_other=other.node_positions.size()+1+last_seen_other;
    }

    auto this_nodes=as_signed_nodes();
    auto other_fw_nodes=other_fw.as_signed_nodes();
    std::vector<sgNodeID_t> merged_nodes;

    //Pick largest left
    if (this_left>=other_left){
        std::copy(this_nodes.begin(),this_nodes.begin()+first_seen_at-1,std::back_inserter(merged_nodes));
    }
    else {
        std::copy(other_fw_nodes.begin(),other_fw_nodes.begin()+first_seen_other-1,std::back_inserter(merged_nodes));
    }

    //merge middle section
    for (auto this_index=first_seen_at-1,other_index=first_seen_other-1; this_index<last_seen_at and other_index<last_seen_other;){
        if (this_nodes[this_index]==other_fw_nodes[other_index]){ //same node on both, add and increment
            merged_nodes.emplace_back(this_nodes[this_index]);
            ++this_index;
            ++other_index;
        }
        else {
            uint64_t this_count=0,other_count=0;
            uint64_t i1;
            for(i1=this_index; i1<last_seen_at and ( disordered_nodes.count(llabs(this_nodes[i1])) or other_fw.node_positions.find(llabs(this_nodes[i1]))==other_fw.node_positions.end()); ++i1) {
                if (disordered_nodes.count(llabs(this_nodes[i1]))==0) ++this_count;
            }
            uint64_t i2;
            for( i2=other_index; i2<last_seen_other and (disordered_nodes.count(llabs(other_fw_nodes[i2])) or node_positions.find(llabs(other_fw_nodes[i2]))==node_positions.end()); ++i2) {
                if (disordered_nodes.count(llabs(other_fw_nodes[i2]))==0) ++other_count;
            }
            if (this_count==0 and other_count==0) { //this is just a single disordered node in both (probably same node reversed?), skip in both
                ++this_index;
                ++other_index;
            }
            else if (this_count>=other_count) {
                for (auto i=this_index;i<i1;++i) {
                    if (disordered_nodes.count(llabs(this_nodes[i]))==0) merged_nodes.emplace_back(this_nodes[i]);
                }
                this_index=i1;
            }
            else {
                for (auto i=other_index;i<i2;++i) {
                    if (disordered_nodes.count(llabs(other_fw_nodes[i]))==0) merged_nodes.emplace_back(other_fw_nodes[i]);
                }
                other_index=i2;
            }



        }
    }

    //Pick largest right
    if (this_right>=other_right){
        std::copy(this_nodes.begin()+last_seen_at,this_nodes.end(),std::back_inserter(merged_nodes));
    }
    else {
        std::copy(other_fw_nodes.begin()+last_seen_other,other_fw_nodes.end(),std::back_inserter(merged_nodes));
    }

    return LocalOrder(merged_nodes);
}

LocalOrder LocalOrder::reverse() const {
    auto nodes=as_signed_nodes();
    std::vector<sgNodeID_t> nodes_rev;
    nodes_rev.reserve(nodes.size());
    for (auto it=nodes.crbegin();it!=nodes.crend();++it) nodes_rev.emplace_back(-*it);
    return LocalOrder(nodes_rev);
}

int64_t LocalOrder::thread_order_crosses(const std::vector<NodePosition> &thread) const {
    sgNodeID_t last_node=0;
    int64_t last_pos=0;
    int64_t crosses=0;
    for (auto &np:thread){
        auto nop=get_node_position(np.node);
        if (nop==0) continue;
        if ( last_node!=0 and ((last_pos>0 != nop>0) or (last_pos>nop))) ++crosses;
        last_node=np.node;
        last_pos=nop;
    }
    return crosses;

}

float HappySorter::node_happiness(sgNodeID_t nid, bool prev, bool next,int min_threads) const {
    if (min_threads==-1) min_threads=min_node_threads;
    uint64_t happy=0,total=0;
    if (nid==0 or rtg.sdg.links.size()-1<llabs(nid)) return 0;
    if (not prev and not next) return 0;
    else if (prev and not next) {
        for (auto &tid:rtg.node_threads(nid,true)){
            if (threads.count(tid)) {
                ++happy;
                ++total;
            }
            else if (rtg.nodes_before_in_thread(tid,nid)>=10){
                ++total;

            }
        }
    }
    else if (not prev and next) {
        for (auto &tid:rtg.node_threads(nid,true)){
            if (threads.count(tid)) {
                ++happy;
                ++total;
            }
            else if (rtg.nodes_after_in_thread(tid,nid)>=10){
                ++total;
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

void HappySorter::recruit_all_happy_threads_q(int min_nodes, int max_span) {
    std::unordered_map<int64_t,int> thread_included_nodes;
    for (const auto &np: order.node_positions) {
        const auto &nid = np.second > 0 ? np.first : -np.first;
        for (auto &tid:rtg.node_threads(nid, true)) ++thread_included_nodes[tid];
    }

    for (auto &tin:thread_included_nodes) {
        if (tin.second>=min_nodes and (thread_included_nodes.count(-tin.first)==0 or thread_included_nodes[-tin.first]<min_nodes) and threads.count(tin.first)==0 and threads.count(-tin.first)==0 and thread_happiness_q(tin.first,min_nodes,max_span)) {
            threads.insert(tin.first);
            fw_open_threads.insert(tin.first);
            bw_open_threads.insert(tin.first);
        }
    }
}

bool HappySorter::thread_happiness_q(int64_t tid, int min_nodes, int max_span) const {
    std::deque<int64_t> hits_pos; //this is overkill, but probably not the performance limiting structure here
    int64_t p=0;
    for (auto np:rtg.get_thread(tid)) {
        ++p;
        if (order.get_node_position(np.node)>0) {
            if (hits_pos.size()==min_nodes) {
                if (p - hits_pos.front() < max_span) return true;
                hits_pos.pop_front();
            }
            hits_pos.emplace_back(p);
        }
    }
    return false;
}

std::unordered_set<sgNodeID_t> HappySorter::find_fw_candidates(float min_happiness, int min_threads, int end_size) const {
    if (min_threads==-1) min_threads=min_node_threads;
    if (min_happiness==-1) min_happiness=min_node_happiness;
    if (end_size==-1) end_size=order_end_size;
    std::unordered_map<sgNodeID_t,int> node_links;
    std::unordered_set<sgNodeID_t> candidates;

    bool debug=false;
    int64_t end_idx=order.size()-end_size;
    if (end_idx<0) end_idx=0;
    if (debug) std::cout<<"looking for nodes with links on "<<fw_open_threads.size()<<" fw_open_threads"<<std::endl;
    if (debug) std::cout<<"end_idx = "<<end_idx<<std::endl;
    for (auto tid:fw_open_threads){
        auto tnps=rtg.get_thread(tid);
        int64_t last_inside_node=0;
        for (;last_inside_node<tnps.size() and order.get_node_position(tnps[last_inside_node].node)>end_idx;++last_inside_node);
        for (auto i=last_inside_node;i<tnps.size();++i) ++node_links[tnps[i].node];
    }
    if (debug) std::cout<<" ... "<<node_links.size()<<" nodes had links"<<std::endl;

    for (auto &nl:node_links) if (nl.second>=min_threads and node_happiness(nl.first,true,false,min_threads)>=min_happiness) candidates.insert(nl.first);
    if (debug) std::cout<<" ... "<<candidates.size()<<" nodes had >="<<min_threads<<" links and >="<<min_happiness<<" happiness"<<std::endl;

    //remove candidates that appear in both directions (should be rare?), also remove candidates located before end_size
    std::set<sgNodeID_t> to_delete;
    for (auto c:candidates) {
        if ( candidates.count(-c)) to_delete.insert(llabs(c));
        else {
            auto p=order.get_node_position(c);
            if (p<0 or (p>0 and p<=order.size()-end_size)) to_delete.insert(c);
        }
    }
    if (debug) std::cout<<" ... "<<to_delete.size()<<" nodes are already present, or candidates in both directions"<<std::endl;
    for (auto cd:to_delete) {
        candidates.erase(cd);
        candidates.erase(-cd);
    }
    if (debug) std::cout<<" ... "<<candidates.size()<<" fw candidates found"<<std::endl;
    return candidates;
}


std::unordered_set<sgNodeID_t> HappySorter::find_bw_candidates(float min_happiness, int min_threads, int end_size) const {
    if (min_threads==-1) min_threads=min_node_threads;
    if (min_happiness==-1) min_happiness=min_node_happiness;
    if (end_size==-1) end_size=order_end_size;
    std::unordered_map<sgNodeID_t,int> node_links;
    std::unordered_set<sgNodeID_t> candidates;
    bool debug=false;
    if (debug) std::cout<<"looking for nodes with links on "<<bw_open_threads.size()<<" bw_open_threads"<<std::endl;
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
//    if (debug) std::cout<<" ... "<<node_links.size()<<" nodes had links"<<std::endl;

    for (auto &nl:node_links) if (nl.second>=min_threads and node_happiness(nl.first,false,true,min_threads)>=min_happiness) candidates.insert(nl.first);
//    if (debug) std::cout<<" ... "<<candidates.size()<<" nodes had >="<<min_threads<<" links and >="<<min_happiness<<" happiness"<<std::endl;
    //remove candidates that appear in both directions (should be rare?), also remove candidates located before end_size
    std::set<sgNodeID_t> to_delete;
    for (auto c:candidates) {
        if ( candidates.count(-c)) to_delete.insert(llabs(c));
        else {
            auto p=order.get_node_position(c);
            if (p<0 or (p>0 and p>=end_size)) to_delete.insert(c);
        }
    }
//    if (debug) std::cout<<" ... "<<to_delete.size()<<" nodes are already present, or candidates in both directions"<<std::endl;
    for (auto cd:to_delete) {
        candidates.erase(cd);
        candidates.erase(-cd);
    }
//    if (debug) std::cout<<" ... "<<candidates.size()<<" bw candidates found"<<std::endl;
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

    //NEW IMPLEMENTATION, to avoid a single hit in the middle of the order to close the threads
//    std::cout<<"thread closing starting, open_fw: "<<fw_open_threads.size()<<", open_bw:"<<bw_open_threads.size()<<std::endl;
    auto nl=order.as_signed_nodes();
    if (nl.size()<order_end) return;
    std::set<sgNodeID_t> bw_nodes,fw_nodes;
    for(auto i=0;i<order_end;++i) bw_nodes.insert(nl[i]);
//    std::cout<<bw_nodes.size()<<" bw_nodes"<<std::endl;
    for (auto tid:bw_open_threads) {
        bool closed=true;
        bool first=true;
        for (auto &p:rtg.get_thread(tid)) {
            if (bw_nodes.count(p.node)) {
                if (not first) closed=false;
                break;
            }
            first=false;
        }
        if (closed) to_close.insert(tid);
    }
//    std::cout<<to_close.size()<<" threads to close after evaluating bw"<<std::endl;
    for(auto i=nl.size()-order_end;i<nl.size();++i) fw_nodes.insert(-nl[i]);
//    std::cout<<fw_nodes.size()<<" fw_nodes"<<std::endl;
    for (auto tid:to_close) bw_open_threads.erase(tid);
    to_close.clear();
    for (auto tid:fw_open_threads) {
        bool closed=true;
        bool first=true;
        for (auto &p:rtg.get_thread(-tid)) { //going in reverse to check last nodes first
            if (fw_nodes.count(p.node)) {
                if (not first) closed=false;
                break;
            }
            first=false;
        }
        if (closed) to_close.insert(tid);
    }
//    std::cout<<to_close.size()<<" threads to close after evaluating fw"<<std::endl;
    for (auto tid:to_close) fw_open_threads.erase(tid);
//    std::cout<<"thread closing finished, open_fw: "<<fw_open_threads.size()<<", open_bw:"<<bw_open_threads.size()<<std::endl;
}

std::map<int64_t, std::vector<std::pair<int64_t, sgNodeID_t>>> HappySorter::make_thread_nodepositions(const std::unordered_set<sgNodeID_t> & nodes, std::set<int64_t> tids) const{
    std::map<int64_t, std::vector<std::pair<int64_t, sgNodeID_t>>> thread_node_positions;

    //make the list of threads
    if (tids.empty()) {
        for (auto nid:nodes) {
            auto ntids = rtg.node_threads(nid, true);
            for (auto tid:ntids) tids.insert(tid);
        }
    }

    //TODO: count all fw/bw in node for each thread and discard those where nodes appear in more than one direciton.
    for (auto &tid:tids){
        auto &tnp=thread_node_positions[tid];
        auto s_nv=rtg.thread_start_nodeview(tid);
        while (order.node_coordinates.count(s_nv.node_id())==0 and nodes.count(s_nv.node_id())==0) {
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
            if (order.node_coordinates.count(ln.node().node_id()) or nodes.count(ln.node().node_id())) {
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


int64_t HappySorter::hs_place_node(const std::unordered_map<sgNodeID_t, int64_t> &node_positions, const std::map<sgNodeID_t, std::vector<std::pair<sgNodeID_t,int64_t>>> & node_distances, sgNodeID_t nid) const {
    if (node_distances.count(nid)==0) return INT64_MIN; //unplaced
    std::vector<int64_t> positions;
    for (const auto & nd:node_distances.at(nid)) {
        if (node_positions.count(nd.first))
            positions.emplace_back(node_positions.at(nd.first)+nd.second);
    }
    if (positions.size()<3) return INT64_MIN;
    std::sort(positions.begin(),positions.end());
//    if (positions[positions.size()-2]-positions[1]>5000) {
//        std::cout<<"position dispersion is too high"<<std::endl;
//        for (auto p:positions) std::cout<<" "<<p;
//        std::cout<<std::endl;
//        //return INT64_MIN;
//    }
    return positions[positions.size()/2]; //Poor man's median hits back
};

void HappySorter::hs_update_npcomplete(std::map<sgNodeID_t, std::pair<bool,bool>> &np_complete,const std::unordered_map<sgNodeID_t, int64_t> &node_positions, const std::map<sgNodeID_t, std::vector<std::pair<sgNodeID_t,int64_t>>> & node_distances, const std::unordered_set<sgNodeID_t> &to_place) const {
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

sgNodeID_t HappySorter::hs_most_connected_node(const std::unordered_map<sgNodeID_t, int64_t> &node_positions, const std::map<sgNodeID_t, std::vector<std::pair<sgNodeID_t,int64_t>>> & node_distances, const std::unordered_set<sgNodeID_t> &to_place) const{
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

std::map<sgNodeID_t, std::vector<std::pair<sgNodeID_t,int64_t>>> HappySorter::hs_tnp_to_distances (const std::map<int64_t, std::vector<std::pair<int64_t, sgNodeID_t>>> &thread_nodepositions,const std::unordered_set<sgNodeID_t> &nodeset) const {
    std::map<sgNodeID_t, std::vector<std::pair<sgNodeID_t,int64_t>>> distances;

    for (const auto &tnp:thread_nodepositions) {
        sgNodeID_t last_node=0;
        int64_t last_pos=0;
        for (auto p:tnp.second){
            if (last_node!=0) {
                if (nodeset.count(p.second)) distances[p.second].emplace_back(last_node,std::max<int64_t>(1,p.first-last_pos));//ensure no two nodes are placed "at the same bp"
                if (nodeset.count(last_node)) distances[last_node].emplace_back(p.second,std::min<int64_t>(-1,last_pos-p.first));//ensure no two nodes are placed "at the same bp"
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
    for (auto &pn:order.node_coordinates) full_nodeset.insert(pn.first);
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
    std::unordered_map<sgNodeID_t, int64_t> node_positions=order.node_coordinates;

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
                        if (verbose) std::cout<<"placed "<<nid<<" @ "<<p<<std::endl;
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
                if (verbose) std::cout<<"placed "<<nid<<" @ "<<p<<" with partial distances"<<std::endl;
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

/**
 * def node_pos(nid,nnp,dists):
    lp=[]
    lpp=[]
    rp=[]
    rp_valid=True
    for a,d in dists[nid]:
        if d>0:
            if a not in nnp:
                #print("nid %d not placed: distance to absent %d is %d" % (nid,a,d))
                return None
            lp.append(nnp[a]+d)
            lpp.append(nnp[a])
        elif a in nnp:
            rp.append(nnp[a]+d)
        else:
            rp_valid=False
    if lp:
        print("placed by lp: ",lp)
        #discard any positions before a node with >1 links to it
        mp=max([x for x,c in Counter(lpp).items() if c>=1])
        lp=[x for x in lp if x>mp]
        return lp[len(lp)//2]
    if not rp_valid:
        return None
    print("placed by rp")
    return rp[len(rp)//2]
 */

int64_t node_pos_ltr(sgNodeID_t nid, const std::unordered_map<sgNodeID_t, int64_t> &nnp, const std::map<sgNodeID_t, std::vector<std::pair<sgNodeID_t,int64_t>>> & dists){
    std::vector<int64_t> lp,rp;
    std::map<int64_t,int64_t> lpp;
    bool rp_valid = true;
    for (auto &ad:dists.at(nid)){
        auto &a=ad.first;auto &d=ad.second;
        if (d>0){
            if (nnp.count(a)==0) return INT64_MIN;
            lp.emplace_back(nnp.at(a)+d);
            ++lpp[nnp.at(a)]; //increase the previous position count
        }
        else if (nnp.count(a)) {
            rp.emplace_back(nnp.at(a)+d);
        }
        else {
            rp_valid = false;
        }
    }
    //placing by previous nodes
    if (not lp.empty()) {
        //discard any proposed positions before the last previous node with >1 links
        int64_t last_pp=INT64_MIN;
        for (auto &pp:lpp) if (pp.second>1 and pp.first>last_pp) last_pp=pp.first;
        std::vector<int64_t> flp;
        for (auto &p:lp) if (p>last_pp) flp.emplace_back(p);
        //poor man's median
        std::sort(flp.begin(),flp.end());
        return flp[flp.size()/2];
    }
    if (not rp_valid) return INT64_MIN;
    std::sort(rp.begin(),rp.end());
    return rp[rp.size()/2];
}
std::vector<std::pair<sgNodeID_t, int64_t>> HappySorter::place_nodes_ltr(const std::unordered_set<sgNodeID_t> &nodes, sgNodeID_t first_node, bool verbose) const {

    if (verbose) std::cout<<"Place nodes LTR starting"<<std::endl;

    //STEP 1: populate node distances for each node
    auto node_distances=hs_tnp_to_distances(make_thread_nodepositions(nodes),nodes);
    if (verbose) std::cout<<"Created distances for "<<node_distances.size()<<" nodes"<<std::endl;

    //STEP 2: place each node in the median of its placed distances.
    std::unordered_map<sgNodeID_t, int64_t> node_positions;
    node_positions[first_node]=0;
    if (verbose) std::cout<<"Placed "<<first_node<<"@ 0"<<std::endl;

    auto to_place=nodes;
    to_place.erase(first_node);


    if (verbose) std::cout<<"Entering placing loop"<<std::endl;
    sgNodeID_t node_to_place;
    int64_t position;
    while (not to_place.empty()){
        position=INT64_MAX;
        for (auto &nid:to_place){
            auto p=node_pos_ltr(nid,node_positions,node_distances);
            if (p!=INT64_MIN and p<position) {
                node_to_place=nid;
                position=p;
            }
        }
        if (position==INT64_MAX) {
            if (verbose) std::cout<<"Can't place any more nodes, aborting"<<std::endl;
            //return {};
            break;
        };//we failed to place any more nodes -> we failed!
        if (verbose) std::cout<<"Placed "<<node_to_place<<"@ "<<position<<std::endl;
        node_positions[node_to_place]=position;
        to_place.erase(node_to_place);
        if (verbose) std::cout<<to_place.size()<<" nodes left to place "<<position<<std::endl;
    }

    //Transform into sorted vector of placed nodes.
    std::vector<std::pair<sgNodeID_t,int64_t>> placements;
    placements.reserve(node_positions.size());
    for (auto &np:node_positions) placements.emplace_back(np.first,np.second);
    std::sort(placements.begin(),placements.end(),[](auto &a,auto &b){return a.second<b.second;});
    return placements;
}

bool HappySorter::add_placed_nodes(const std::vector<std::pair<sgNodeID_t, int64_t>> &placed_nodes,
                                   bool update_current) {
    auto old_node_count=order.node_positions.size();
    for (auto pn:placed_nodes) {
        if (update_current or order.node_coordinates.count(pn.first)==0) order.node_coordinates[pn.first]=pn.second;
    }
    //std::vector<std::pair<sgNodeID_t, int64_t>> all_nodes(order.node_coordinates.begin(),order.node_coordinates.end());
    std::vector<std::pair<sgNodeID_t, int64_t>> all_nodes;
    for (auto &nc:order.node_coordinates) all_nodes.emplace_back(nc.first,nc.second);
    std::sort(all_nodes.begin(),all_nodes.end(),[](auto &a,auto &b){return a.second<b.second;});
    order.node_positions.clear();
    for (auto i=0;i<all_nodes.size();++i) order.node_positions[llabs(all_nodes[i].first)]=all_nodes[i].first>0? i+1:-i-1;
    return order.node_positions.size()>old_node_count;
}

bool HappySorter::q_grow_loop(int min_threads, float min_happiness, int p, int q, int64_t steps, bool verbose) {
    if (min_threads==-1) min_threads=min_node_threads;
    if (min_happiness==-1) min_happiness=min_node_happiness;
    bool grown=false;
    int64_t step=0;
    while (++step<=steps) {
        /********/

        auto old_size=order.node_coordinates.size();
        std::unordered_set<sgNodeID_t> candidates;
//        sdglib::OutputLog()<<"adding FW candidates"<<std::endl;

        auto pn=place_nodes(find_fw_candidates(min_happiness, min_threads), false);
        { //XXX: temporary fix that discards nodes that are placed before nodes included in the order already
            int64_t first_included = INT64_MAX;
            std::vector<sgNodeID_t> to_delete;
            for (auto &p:pn) if (order.node_coordinates.count(p.first) and first_included>p.second)first_included=p.second;
            for (auto &p:pn) if (order.node_coordinates.count(p.first) and first_included>p.second) to_delete.push_back(p.first);
            if (not to_delete.empty()) {
                std::cout<<"WARNING: removing "<<to_delete.size()<<" fw placed nodes placed before the end nodes start"<<std::endl;
            }
        }
        add_placed_nodes(pn);
//        sdglib::OutputLog()<<"adding BW candidates"<<std::endl;
//        add_placed_nodes(place_nodes(find_bw_candidates(min_happiness, min_threads), false));
        pn=place_nodes(find_bw_candidates(min_happiness, min_threads), false);
        { //XXX: temporary fix that discards nodes that are placed after nodes included in the order already
            int64_t first_included = INT64_MIN;
            std::vector<sgNodeID_t> to_delete;
            for (auto &p:pn) if (order.node_coordinates.count(p.first) and first_included<p.second)first_included=p.second;
            for (auto &p:pn) if (order.node_coordinates.count(p.first) and first_included<p.second) to_delete.push_back(p.first);
            if (not to_delete.empty()) {
                std::cout<<"WARNING: removing "<<to_delete.size()<<" bw placed nodes placed after the end nodes start"<<std::endl;
            }
        }
        add_placed_nodes(pn);
        auto ends_to_fill = 2*(order.node_coordinates.size()-old_size);
        if (ends_to_fill<200) ends_to_fill=200;
//        sdglib::OutputLog()<<"adding internal candidates"<<std::endl;
        add_placed_nodes(place_nodes(find_internal_candidates(min_happiness, min_threads,0,ends_to_fill), false));
        add_placed_nodes(place_nodes(find_internal_candidates(min_happiness, min_threads,order.node_coordinates.size()-ends_to_fill,order.node_coordinates.size()), false));
//        sdglib::OutputLog()<<"updating positions"<<std::endl;
        update_positions(0,ends_to_fill);
        update_positions(order.node_coordinates.size()-ends_to_fill,order.node_coordinates.size());

//        sdglib::OutputLog()<<"recruiting threads"<<std::endl;
        recruit_all_happy_threads_q(p,q);
//        sdglib::OutputLog()<<"closing threads"<<std::endl;
        close_internal_threads();

        grown=order.node_coordinates.size()>old_size;
        /********/

        if (verbose and (step%50==0 or step==steps or not grown)) {
            std::pair<sgNodeID_t,int64_t> cmax={0,0},cmin={0,0};
            for (const auto &nc:order.node_coordinates) {
                if (nc.second<cmin.second) cmin=nc;
                if (nc.second>cmax.second) cmax=nc;
            }
            sdglib::OutputLog()<<"After "<<step<<" grow rounds: "<<order.node_coordinates.size()<<" anchors span "<<cmax.second-cmin.second<<" bp, "<<cmin.first<<" @ "<<cmin.second<<" : "<<cmax.first<<" @ "<<cmax.second<<std::endl;
        }
        if (not grown) break;
    }
    return grown;
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

bool HappySorter::update_positions(int64_t first, int64_t last) {
//    std::cout<<std::endl<<"HS::update_positions starting"<<std::endl;
    if(last==-1 or last>order.size()-1) last=order.size()-1;
    if (first>order.size() or first>last) return false;

    //First, create nodepositions and distances, but only on threads from the nodes to update.
    auto onodes=order.as_signed_nodes();

    std::unordered_set<sgNodeID_t> nodes(onodes.cbegin()+first,onodes.cbegin()+last+1);
    std::unordered_set<sgNodeID_t> all_nodes(onodes.cbegin(),onodes.cend());

    std::set<int64_t> tids;
    for (auto nid:nodes) {
        auto ntids = rtg.node_threads(nid, true);
        for (auto tid:ntids) tids.insert(tid);
    }

    auto tnp=make_thread_nodepositions(all_nodes,tids);
    auto nd=hs_tnp_to_distances(tnp,all_nodes);

    //sort the nodes to update by their distance to the origin
    std::vector<sgNodeID_t> nodes_to_update(onodes.cbegin()+first,onodes.cbegin()+last+1);
    std::sort(nodes_to_update.begin(),nodes_to_update.end(),[this](auto a, auto b){return llabs(order.node_coordinates[a])<llabs(order.node_coordinates[b]);});
//    std::cout<<std::endl<<nodes_to_update.size()<<" nodes to update from "<<all_nodes.size()<<" nodes in the order"<<std::endl;

    bool changed=true;

    bool ever_changed=false;
    for (auto steps=10; steps and changed;--steps) {
//        std::cout<<std::endl<<"HS::update_positions starting reposition round"<<std::endl;
        std::vector<sgNodeID_t> to_place;
        changed=false;
        //update node positions for any nodes with >3 links, if less, skip them to be fully placed later
        std::vector<int64_t> nps;
        for (auto &nid:nodes_to_update){
            //if (order.node_coordinates[nid]==0) continue;
            nps.clear();
            for (auto &d:nd[nid]) {
                if ((d.second>0)==(order.node_coordinates[nid]>0)) nps.emplace_back(order.node_coordinates[d.first]+d.second);
            }
            if (nps.empty()) {
                to_place.emplace_back(nid);
                continue;
            }
//            if (nps.empty()) continue;
            std::sort(nps.begin(),nps.end());
            auto new_p=nps[nps.size()/2];
            if (new_p*order.node_coordinates[nid]>0 and new_p != order.node_coordinates[nid]) {
                changed= true;
//                std::cout<<"moving "<<nid<<" from "<<order.node_coordinates[nid]<<" to "<<new_p<<std::endl;
                order.node_coordinates[nid]=new_p;
            }
        }
        //place nodes with <=3 links
        for (auto nid:to_place) {
            auto new_p=hs_place_node(order.node_coordinates,nd,nid);
            if (new_p != order.node_coordinates[nid]) {
                changed= true;
                order.node_coordinates[nid]=new_p;
            }
        }
        if (changed) ever_changed=true;
//        if (onodes.size()!=order.node_coordinates.size()) std::cout<<"coordinates ("<<order.node_coordinates.size()<<") and signed nodes ("<<onodes.size()<<") have different sizes! "<<std::endl;
    }

    //re-place all nodes
    for (auto steps=10; steps and changed;--steps) {
        changed = false;
        for (auto nid:nodes_to_update) {
            auto new_p = hs_place_node(order.node_coordinates, nd, nid);
            if (new_p != order.node_coordinates[nid]) {
                changed = true;
                order.node_coordinates[nid] = new_p;
            }
        }
        if (changed) ever_changed = true;
    }
//    std::cout<<std::endl<<"HS::update_positions updating order"<<std::endl;
    std::vector<std::pair<sgNodeID_t,int64_t>> ordered_nodecoords(order.node_coordinates.begin(),order.node_coordinates.end());
    std::sort(ordered_nodecoords.begin(),ordered_nodecoords.end(),[](auto a, auto b){return a.second<b.second;});
//    std::cout<<ordered_nodecoords.size()<<" nodes by coordinate, and "<<order.node_positions.size()<<" node positions"<<std::endl;
//    if (onodes.size()!=ordered_nodecoords.size()) std::cout<<"nodes by coordinate ("<<ordered_nodecoords.size()<<") and signed nodes ("<<onodes.size()<<") have different sizes! "<<std::endl;
    for (auto i=0;i<ordered_nodecoords.size();++i) {
        const auto &np=ordered_nodecoords[i];
        if (np.first>0 and order.node_positions[np.first]!=i+1) {
//            std::cout<<"moving "<<np.first<<" from "<<order.node_positions[np.first]<<" to "<<i+1<<std::endl;
            order.node_positions[np.first]=i+1;
        }
        else if (np.first<0 and order.node_positions[-np.first]!=-i-1){
//            std::cout<<"moving "<<-np.first<<" from "<<order.node_positions[-np.first]<<" to "<<-i-1<<std::endl;
            order.node_positions[-np.first]=-i-1;
        }
    }
//    std::cout<<"after updating order: "<<ordered_nodecoords.size()<<" nodes by coordinate, and "<<order.node_positions.size()<<" node positions"<<std::endl;
//    std::cout<<std::endl<<"HS::update_positions finished"<<std::endl;
    return ever_changed;
}

void HappySorter::run_from_nodelist(std::vector<sgNodeID_t> nodes, int min_threads, float min_happiness, int p, int q,
                                    int64_t steps, bool verbose) {
    if (min_threads==-1) min_threads=min_node_threads;
    if (min_happiness==-1) min_happiness=min_node_happiness;
    std::vector<int> thread_counts;
    std::vector<sgNodeID_t> usable_nodes;
    for (auto &n:nodes) {
        auto ntc=rtg.node_threads(n).size();
        if (ntc>=min_threads*4) {
            usable_nodes.push_back(n);
            thread_counts.push_back(ntc);
        }
    }
    std::sort(thread_counts.begin(),thread_counts.end());
    auto median_tc=thread_counts[thread_counts.size()/2];
    for (auto &n:usable_nodes){
        if (rtg.node_threads(n).size()>=median_tc) {
            //Force thread recruiting from the current nodes
            if (verbose) sdglib::OutputLog()<<"Starting from node "<<n<<std::endl;
            threads.clear();
            fw_open_threads.clear();
            bw_open_threads.clear();
            order=LocalOrder(usable_nodes);
            if (verbose) sdglib::OutputLog()<<"Local order filled temporally with all "<<order.size()<<" usable nodes"<<std::endl;
            recruit_all_happy_threads_q(p,q);
            if (verbose) sdglib::OutputLog()<<"Threads recruited: "<<threads.size()<<std::endl;
            order=LocalOrder();
            if (verbose) sdglib::OutputLog()<<"Starting from node "<<n<<std::endl;
            add_placed_nodes({{n,0}});

            q_grow_loop(min_threads,min_happiness,p,q,steps,verbose);
            if (order.size()>usable_nodes.size()) break; //keeps trying until it grows the order, or runs out of nodes
        }
    }
}


void HappySorterRunner::run_from_dg_lines(DistanceGraph dg, int min_line_nodes, int p, int q) {
    sdglib::OutputLog()<<"Starting HappySorterRunner run from rtg lines..."<<std::endl;
    //sdglib::OutputLog()<<"Getting all nodeviews..."<<std::endl;
    auto lines=dg.get_all_lines(min_line_nodes);
    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    shuffle (lines.begin(), lines.end(), std::default_random_engine(seed));
    sdglib::OutputLog()<<lines.size()<<" lines shuffled..."<<std::endl;
#pragma omp parallel for schedule(dynamic,1)
    for (auto lineit=lines.cbegin();lineit<lines.cend();++lineit){
        auto &line=*lineit;
        try {
            auto hs = HappySorter(rtg, min_thread_happiness, min_node_happiness, min_thread_nodes, min_node_threads, order_end_size);
            hs.run_from_nodelist(line,min_node_threads,min_node_happiness,p,q);
#pragma omp critical(orders)
            {
                orders[orders.size()]=hs.order;
                if (orders.size()%10==0) sdglib::OutputLog()<<orders.size()<<" orders created"<<std::endl;
            }

        }
        catch (std::exception &e) {
            std::cout<<"Exception caught when trying to process line with node "<<line[0]<<": "<< typeid(e).name() <<": "<<e.what()<<", skipping"<<std::endl;
        }

    }
    sdglib::OutputLog()<<"HappySorterRunner run finished!!!"<<std::endl;
}