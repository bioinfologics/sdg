//
// Created by Bernardo Clavijo (EI) on 06/08/2021.
//

#include "RTGPartition.hpp"

template <typename T>
std::vector<T> find_all_p_in_q(int p, int q, const std::vector<T>& input, bool discard_zero=false){
    std::vector<T> v;
    std::unordered_map<T,int64_t> votes;
    for (auto &x:input) if (not discard_zero or x!=0) ++votes[x];
    std::vector<int64_t> pos;
    for (auto xc:votes) if (xc.second>=p) {
        for (auto i=0;i<input.size();++i) {
            auto x=input[i];
            if (x == xc.first) {
                pos.emplace_back(i);
                if (pos.size()>=p and pos.back()-*(pos.end()-p)<q) {
                    v.emplace_back(x);
                    break;
                }
            }
        }
        pos.clear();
    }
    return v;
};

//TODO: this needs parallelisation. fill things in local maps of vectors, then swap them into the real maps.
RTGPartition::RTGPartition(const ReadThreadsGraph &rtg, int min_node_threads, float node_min_percentage, int thread_p,
                             int thread_q, int min_thread_nodes): rtg(rtg), min_node_threads(min_node_threads),
                             node_min_percentage(node_min_percentage), thread_p(thread_p), thread_q(thread_q), min_thread_nodes(min_thread_nodes) {
    if (min_thread_nodes==-1) this->min_thread_nodes=thread_p;
    sdglib::OutputLog()<<"analysing nodeviews"<<std::endl;
    for (auto &nv:rtg.get_all_nodeviews(false,false)) {
        auto nid=nv.node_id();
        auto nts=rtg.node_threads(nid);
        node_threads[nid].reserve(nts.size());
        node_class[nid]=0;
        for (auto &tid : nts) {
            if (rtg.thread_info.at(tid).link_count<=this->min_thread_nodes-1) continue;
            node_threads[nid].emplace_back(tid);
            thread_class[tid]=0;
        }
        if (node_threads[nid].empty()) node_threads.erase(nid);
    }
    sdglib::OutputLog()<<"analysing threads"<<std::endl;
    for (auto & tc:thread_class){
        auto t=rtg.get_thread_nodes(tc.first);
        thread_nodes[tc.first].reserve(t.size());
        for (auto &nid:t) thread_nodes[tc.first].emplace_back(llabs(nid));
    }
    sdglib::OutputLog()<<"done with "<<node_class.size()<<" nodes in "<<thread_class.size()<<" threads"<<std::endl;
}


int64_t RTGPartition::get_node_class(sgNodeID_t nid) {
    if (node_class.count(nid)==0) return -1;
    return node_class[nid];
}
int64_t RTGPartition::compute_node_class(sgNodeID_t nid) {
    std::map<int64_t, int64_t> class_votes;
    int64_t highest_votes=0;
    int64_t highest_votes_class=0;
    int64_t v;
    int64_t c;
    auto &threads=node_threads[nid];
    for (auto tid:threads) {
        c= thread_class[tid];
        if (c==0) continue;
        v = ++class_votes[c];
        if (highest_votes < v) {
            highest_votes = v;
            highest_votes_class = c;
        }
    }
    if (highest_votes>=min_node_threads and highest_votes>=node_min_percentage*threads.size()) {
        for (auto &cv:class_votes) {
            if (cv.first!=0 and cv.first!=highest_votes_class and cv.second/4>=highest_votes/*node_min_percentage*threads.size()*/) return 0; //two classes met the number of votes
        }
        return highest_votes_class;
    }
    return 0;
}

bool RTGPartition::switch_node_class(sgNodeID_t nid, int64_t c) {
    if (node_class.count(nid)==0 or node_class[nid]==c) return false;
    node_class[nid]=c;
    for (auto &tid:node_threads[nid]) threads_to_evaluate.insert(tid);
    return true;
}

int64_t RTGPartition::get_thread_class(int64_t tid) {
    if (thread_class.count(tid)==0) return -1;
    return thread_class[tid];
}

int64_t RTGPartition::compute_thread_class(int64_t tid, int distance_to_end) {
    std::map<int64_t, int64_t> class_votes;
    int64_t highest_votes=0;
    int64_t highest_votes_class=0;


    auto &nodes=thread_nodes[tid];
    std::vector<int64_t> nc;
    nc.reserve(nodes.size());
    for (auto nid:nodes) nc.emplace_back(node_class[nid]);

    int64_t v;
    for (auto c:nc) {
        if (c==0) continue;
        v = ++class_votes[c];
        if (highest_votes < v) {
            highest_votes = v;
            highest_votes_class = c;
        }
    }
    int64_t winner_class=0;
    if (2*highest_votes>nodes.size()) //condition 1: 50+% nodes in class
        winner_class=highest_votes_class;
    else if (highest_votes>=thread_p) { //condition 2: only one class with p in q
        auto p_in_q=find_all_p_in_q(thread_p,thread_q,nc,true);
        if (p_in_q.size()==1) winner_class=p_in_q[0];
    }
//    if (winner_class!=0 and nc.size()>distance_to_end) {
//        bool start=false;
//        bool end=false;
//        for (auto i=0;i<distance_to_end;++i) {
//            if (nc[i]==winner_class) {
//                start=true;
//                break;
//            }
//        }
//        for (auto i=nc.size()-distance_to_end ; i<nc.size();++i) {
//            if (nc[i]==winner_class) {
//                end=true;
//                break;
//            }
//        }
//        if (not start and not end) winner_class=0;
//    }
    return winner_class;
}

bool RTGPartition::switch_thread_class(int64_t tid, int64_t c) {
    if (thread_class.count(tid)==0 or thread_class[tid]==c) return false;
    thread_class[tid]=c;
    for (auto &nid:thread_nodes[tid]) nodes_to_evaluate.insert(nid);
    return true;
}

std::vector<int64_t> RTGPartition::find_unclassified_threads(int min_nodes, float max_classified_nodes_perc) {
    std::vector<int64_t> tids;
    for (auto tc:thread_class) {
        if (thread_available(tc.first,min_nodes,max_classified_nodes_perc)) tids.emplace_back(tc.first);
    }
    return tids;
}

int64_t RTGPartition::new_class_from_thread(int64_t tid) {
    //find the next-up class
    int64_t next_class=1;
    for (auto &tc:thread_class) if (tc.second>=next_class) next_class=tc.second+1;
    for (auto &nc:node_class) if (nc.second>=next_class) next_class=nc.second+1;
    switch_thread_class(tid,next_class);
    int64_t node_counter=0;
    for (auto nid:thread_nodes[tid]){
        if (node_class[nid]==0) {
            ++node_counter;
            switch_node_class(nid,next_class);
        }
    }
    return next_class;
}

LocalOrder RTGPartition::order_from_class(int64_t cid) {

    std::unordered_set<sgNodeID_t> nodes_to_place,nodes;
    std::unordered_set<int64_t> threads_to_place,threads;
    for (auto nc:node_class) if (nc.second==cid) nodes_to_place.insert(nc.first);
    for (auto tc:thread_class) if (tc.second==cid) threads_to_place.insert(tc.first);

    std::unordered_map<sgNodeID_t,int64_t> node_votes;
    std::unordered_map<int64_t,int64_t> thread_votes;

    std::cout<<"placing a node and its threads"<<std::endl;
    //choose a node, place it in fw orientation, place all its threads
    sgNodeID_t starting_node=0;
    for (auto nid:nodes_to_place) {
        if (node_threads[nid].size()>6) {
            nodes.insert(nid);
            starting_node=nid;
            nodes_to_place.erase(nid);
            for (auto tid:rtg.node_threads(nid,true)) if (threads_to_place.count(llabs(tid))) {
                threads.insert(tid);
                for (auto nid:rtg.get_thread_nodes(tid,true)) if (nodes_to_place.count(llabs(nid))) ++node_votes[nid];
            }
            break;
        }
    }
    uint64_t last_to_place=nodes_to_place.size()+threads_to_place.size()+1;

    while (nodes_to_place.size()+threads_to_place.size()!=last_to_place) {
//        std::cout<<"Starting placement round, "<<nodes_to_place.size()+threads_to_place.size()<<" elements to place, previous count was "<<last_to_place<<std::endl;
        last_to_place=nodes_to_place.size()+threads_to_place.size();
        std::unordered_set<sgNodeID_t> new_nodes;
        std::unordered_set<int64_t> new_threads;
//        std::cout<<"placing nodes"<<std::endl;
        for (auto nid:nodes_to_place) {
            for (auto nnid:{nid,-nid}){
                if (node_votes[nnid]>=2) {
                    new_nodes.insert(nnid);
                    nodes.insert(nnid);
                    for (auto tid:rtg.node_threads(nnid,true)) if (threads_to_place.count(llabs(tid))) ++thread_votes[tid];
                    break;
                }
            }
        }
//        std::cout<<"placing threads"<<std::endl;
        for (auto tid:threads_to_place) {
            for (auto ttid:{tid,-tid}){
                if (thread_votes[ttid]>=2) {
                    new_threads.insert(ttid);
                    threads.insert(ttid);
                    for (auto nid:rtg.get_thread_nodes(ttid,true)) if (nodes_to_place.count(llabs(nid))) ++node_votes[nid];
                    break;
                }
            }
        }
//        std::cout<<"removing "<<new_nodes.size()<<" placed nodes from to_place"<<std::endl;
        for (auto nid:new_nodes) nodes_to_place.erase(llabs(nid));
//        std::cout<<"removing "<<new_threads.size()<<" placed threads from to_place"<<std::endl;
        for (auto tid:new_threads) threads_to_place.erase(llabs(tid));
//        std::cout<<"placement round finished, "<<nodes_to_place.size()<<" nodes, and "<<threads_to_place.size()<<" threads left to place"<<std::endl;
    }

    HappySorter hs(rtg,node_min_percentage,min_node_threads,0); //TODO options for happiness, etc
    //now add all threads
    hs.order.add_placed_nodes({{starting_node,0}});
    hs.threads=threads;
    //now place every node
    //now place a random node as start
    hs.order.add_placed_nodes(hs.place_nodes(nodes,false));
    return LocalOrder(hs.order);
}


void RTGPartition::classify_all_threads(int min_nodes, float max_classified_nodes_perc) {
    auto ut=find_unclassified_threads(min_nodes,max_classified_nodes_perc);
    for (auto i=0;i<ut.size();++i) {
        auto tid=ut[i];
        if (thread_available(tid,min_nodes,max_classified_nodes_perc)) {
            //std::cout<<"starting sorter from thread"<<tid<<std::endl;
            auto new_class=new_class_from_thread(tid);
            int64_t nodes_in_class=0;
            for (const auto & nc:node_class) if (nc.second==new_class) ++nodes_in_class;
            //std::cout<<"Class "<<new_class<<" from sorter from thread "<<tid<<" added with "<<nodes_in_class<<" / "<<hs.order.size()<<" nodes, propagating..."<<std::endl;
            propagate();
            nodes_in_class=0;
            for (const auto & nc:node_class) if (nc.second==new_class) ++nodes_in_class;
            std::cout<<"Class "<<new_class<<" from sorter from thread "<<tid<<" has "<<nodes_in_class<<" nodes after propagation"<<std::endl;
            auto classified_nodes=0;
            for (const auto &nc:node_class) if (nc.second!=0) ++classified_nodes;
            std::cout<<classified_nodes<<" / "<<node_class.size()<<" nodes classified"<<std::endl;
        }
    }
}

uint64_t RTGPartition::propagate(uint64_t steps, bool verbose) {
    uint64_t switched=0;
    while (steps-- and threads_to_evaluate.size()+nodes_to_evaluate.size()) {
        std::unordered_set<int64_t> otte;
        std::swap(otte,threads_to_evaluate);
        for (auto tid:otte) {
            if(switch_thread_class(tid,compute_thread_class(tid))) {
                ++switched; //switch won't do anything if class doesn't change
                if (verbose) std::cout<<"Thread "<<tid<<" switched class to "<<thread_class[tid]<<std::endl;
            }
        }
        std::unordered_set<sgNodeID_t> onte;
        std::swap(onte,nodes_to_evaluate);
        for (auto nid:onte){
            if (switch_node_class(nid,compute_node_class(nid))) {
                ++switched; //switch won't do anything if class doesn't change
                if (verbose) std::cout<<"Node "<<nid<<" switched class to "<<node_class[nid]<<std::endl;
            }
        }
    }
    return switched;
}

uint64_t RTGPartition::thread_propagate(uint64_t steps, float vote_perc, int max_noise, bool verbose) {

    uint64_t switched=0;
    while (steps--){
        std::map<sgNodeID_t, std::map<int64_t,uint64_t>> tcvotes;
        std::vector<std::pair<uint64_t,uint64_t>> to_switch;
        if (verbose) std::cout<<"thread propagation starting, "<<steps+1<<" steps left"<<std::endl;
        //go through every thread that has a class, check its neighbours for p/q condition, if met, add a vote to switch them to this class
        for (auto &tc:thread_class){
            if (tc.second!=0) {
                for (auto &nt:thread_neighbours[tc.first]) {
                    ++tcvotes[nt][tc.second];
                }
            }
        }
        if (verbose) std::cout<<tcvotes.size()<<" threads with votes"<<std::endl;
        //check votes, if a thread has only one class with enough votes and no other class with >1 vote, switch it
        uint64_t new_switched=0;
        for (auto &tcv:tcvotes){
            int64_t winner=0;
            uint64_t threshold=((uint64_t) (vote_perc * 100 )) * thread_neighbours[tcv.first].size() / 100;
            threshold=(2*max_noise>threshold? 2*max_noise:threshold);
            for (auto &cv:tcv.second) {
                if (winner==0 and cv.second>=threshold){
                    winner=cv.first;
                }
                else if (cv.second>max_noise) {
                    winner=0;
                    break;
                }

            }
            if (winner!=thread_class[tcv.first]) {
//                if (verbose) {
//                    std::cout<<"Thread "<<tcv.first<<" to class "<<winner<<", threshold "<<threshold<<" out of "<<thread_neighbours[tcv.first].size()<<" neighbours, votes:";
//                    for (auto &cv:tcv.second) std::cout<<" "<<cv.first<<": "<<cv.second<<",";
//                    std::cout<<std::endl;
//                    to_switch.emplace_back(tcv.first,winner);
//                }
                ++new_switched;
            }
        }
        for (auto tc:to_switch) switch_thread_class(tc.first,tc.second);
        if (verbose) std::cout<<new_switched<<" threads switched"<<std::endl;
        if (new_switched==0) break;//no switches, propagation finished.
        switched+=new_switched;
    }
    return switched;
}



uint64_t RTGPartition::thread_propagate2(uint64_t steps, int min_side_votes, float side_vote_perc,
                                          int min_contain_votes, float contain_vote_perc, bool reclassify, bool verbose) {


    uint64_t switched=0;
    while (steps--){
        //Steps: 1 -> find candidates  / 3 -> find candidates / 4 -> recruit fw/bw threads (check coherency first?)
        std::set<uint64_t> class_neighbours;
        for (auto &tn:thread_neighbours){
            auto c=thread_class[tn.first];
            for (auto &ntid:tn.second) {
                if (thread_class[ntid]!=c) {
                    class_neighbours.insert(tn.first);
                    break;
                }
            }
        }
        if (verbose) std::cout<<"there are "<<class_neighbours.size()<<" class neighbours"<<std::endl;
        // 2 -> recruit threads
        std::vector<std::pair<uint64_t,uint64_t>> to_switch;
        for (auto tid:class_neighbours) {
            if (not reclassify and thread_class[tid]!=0) continue;
            auto &neighs=thread_neighbours[tid];
            auto &neightypes=thread_neighbours_types[tid];
            if (neighs.size()!=neightypes.size()) continue;
            //TODO: add an "included" case, where all threads included in this one vote to place it
            uint64_t including_total=0,fw_total=0,bw_total=0;
            uint64_t including_highest=0,fw_highest=0,bw_highest=0;
            uint64_t including_highest_class=0,fw_highest_class=0,bw_highest_class=0;
            std::map<uint64_t,uint64_t> included_votes,bw_votes,fw_votes;
            for (auto nix=0;nix<neighs.size();++nix){
                auto ntid=neighs[nix];
                auto ntype=neightypes[nix];
                //count included votes
                if (ntype==ThreadOverlapType::t1_in_t2 or ntype==ThreadOverlapType::t1_in_rt2 or ntype==ThreadOverlapType::complete) {
                    ++including_total;
                    auto c= thread_class[ntid];
                    if (c!=0) {
                        auto v = ++included_votes[c];
                        if (including_highest < v) {
                            including_highest = v;
                            including_highest_class = c;
                        }
                    }
                }
                else if (ntype==ThreadOverlapType::t1_t2 or ntype==ThreadOverlapType::t1_rt2) {
                    ++fw_total;
                    auto c= thread_class[ntid];
                    if (c!=0) {
                        auto v = ++fw_votes[c];
                        if (fw_highest < v) {
                            fw_highest = v;
                            fw_highest_class = c;
                        }
                    }
                }
                else if (ntype==ThreadOverlapType::t2_t1 or ntype==ThreadOverlapType::rt2_t1) {
                    ++fw_total;
                    auto c= thread_class[ntid];
                    if (c!=0) {
                        auto v = ++fw_votes[c];
                        if (fw_highest < v) {
                            fw_highest = v;
                            fw_highest_class = c;
                        }
                    }
                }

            }
            //TODO: consider the winners, if clear cut (all winners to same class, where they win), add to to_switch
            int64_t winner=0;
            if (including_highest>=min_contain_votes and including_highest>=including_total*contain_vote_perc)
                winner=including_highest_class;
            if (fw_highest>=min_side_votes and fw_highest>=fw_total*side_vote_perc) {
                if (winner==0 or winner==fw_highest_class) winner=fw_highest_class;
                else winner=-1;
            }
            if (bw_highest>=min_side_votes and bw_highest>=bw_total*side_vote_perc) {
                if (winner==0 or winner==bw_highest_class) winner=bw_highest_class;
                else winner=-1;
            }
            std::cout<<tid<<": IN: "<<including_highest_class<<" ("<<including_highest<<"/"<<including_total<<")  BW: "<<bw_highest_class<<" ("<<bw_highest<<"/"<<bw_total<<")  FW:"<<fw_highest_class<<" ("<<fw_highest<<"/"<<fw_total<<") --> winner: "<<winner<<std::endl;
            if (winner>0 and winner!=thread_class[tid]) {
                std::cout<<" SWITCHED from "<<thread_class[tid]<<" to "<<winner<<"!"<<std::endl;
                to_switch.emplace_back(tid,winner);
            }
        }

        for (auto tc:to_switch) switch_thread_class(tc.first,tc.second);
        if (verbose) std::cout<<to_switch.size()<<" threads switched"<<std::endl;
        if (to_switch.empty()) break;//no switches, propagation finished.
        switched+=to_switch.size();
    }
    return switched;
}

std::map<uint64_t, uint64_t> RTGPartition::thread_class_votes(uint64_t tid, std::set<ThreadOverlapType> ovltypes) {
    std::map<uint64_t, uint64_t> votes;
    auto &tn=thread_neighbours[tid];
    auto &tnt=thread_neighbours_types[tid];
    for (auto i=0;i<tn.size();++i) {
        if (ovltypes.empty() or ovltypes.count(tnt[i])) ++votes[thread_class[tn[i]]];
    }
    return votes;
}

std::unordered_map<std::pair<int64_t, int64_t>, std::vector<int64_t>> RTGPartition::find_class_bridges(int p, int q) {
    std::unordered_map<std::pair<int64_t, int64_t>, std::vector<int64_t>> cb;
    for (auto tn:thread_nodes){
        //if (thread_class[tn.first]!=0) continue;
        auto & nodes=tn.second;
        std::vector<int64_t> nc;
        nc.reserve(nodes.size());
        for (auto nid:nodes) nc.emplace_back(node_class[nid]);
        auto pq=find_all_p_in_q(p,q,nc,true);
        if (pq.size()==2) {
            cb[{std::min(pq[0],pq[1]),std::max(pq[0],pq[1])}].emplace_back(tn.first);
        }
    }
    return cb;
}

void RTGPartition::compute_thread_intersections(int min_threads, int max_threads) {
    thread_intersections.clear();
    for (auto &nt:node_threads){
        if (nt.second.size()<min_threads or nt.second.size()>max_threads) continue;
        for (auto i=0;i<nt.second.size()-1;++i){
            for (auto j=i;j<nt.second.size();++j){
                if (nt.second[i]<nt.second[j]) ++thread_intersections[{nt.second[i],nt.second[j]}];
                else ++thread_intersections[{nt.second[j],nt.second[i]}];
            }
        }
    }
}

void RTGPartition::compute_thread_neighbours(int min_shared) {
    thread_neighbours.clear();
    for (auto &ti:thread_intersections) {
        if (ti.second>=min_shared){
            thread_neighbours[ti.first.first].emplace_back(ti.first.second);
            thread_neighbours[ti.first.second].emplace_back(ti.first.first);
        }
    }
}

void RTGPartition::compute_thread_neighbours_p_q(int p, int q) {
    thread_neighbours.clear();
#pragma omp parallel for
    for (auto b=0;b<thread_intersections.bucket_count();++b) {
        for(auto bi=thread_intersections.begin(b);bi!=thread_intersections.end(b);bi++){
            auto &ti=*bi;
            auto t1 = ti.first.first;
            auto t2 = ti.first.second;
            if (t1 == t2) continue; //self intersection is there only to show a total of analysed nodes in the thread
            if (ti.second >= p and find_all_p_in_q(p, q, thread_shared_detail(t1, t2), true).size() == 1 and
                find_all_p_in_q(p, q, thread_shared_detail(t2, t1), true).size() == 1) {
#pragma omp critical
                {
                    thread_neighbours[t1].emplace_back(t2);
                    thread_neighbours[t2].emplace_back(t1);
                }
            }
        }

    }
}

std::vector<int64_t> RTGPartition::get_thread_neighbours(int64_t tid) const {
    if (thread_neighbours.count(llabs(tid))==0) return {};
    return thread_neighbours.at(llabs(tid));
}

std::vector<ThreadOverlapType> RTGPartition::get_thread_neighbours_types(int64_t tid) const {
    if (thread_neighbours_types.count(llabs(tid))==0) return {};
    return thread_neighbours_types.at(llabs(tid));
}

void RTGPartition::reset(int _min_node_threads, float _node_min_percentage, int _thread_p, int _thread_q) {
    if (_min_node_threads!=-1) min_node_threads=_min_node_threads;
    if (_node_min_percentage!=-1) node_min_percentage=_node_min_percentage;
    if (_thread_p!=-1) thread_p=_thread_p;
    if (_thread_q!=-1) thread_p=_thread_q;
    for (auto &tc:thread_class) tc.second=0;
    for (auto &nc:node_class) nc.second=0;
    nodes_to_evaluate.clear();
    threads_to_evaluate.clear();
    thread_intersections.clear();
}

int64_t RTGPartition::get_thread_intersection(int64_t tid1, int64_t tid2) const {
    std::pair<int64_t, int64_t> p;
    tid1 = llabs(tid1);
    tid2 = llabs(tid2);
    if (tid1<tid2) p={tid1,tid2};
    else p={tid2,tid1};
    if (thread_intersections.count(p)) return thread_intersections.at(p);
    return 0;
}

std::vector<int> RTGPartition::thread_shared_detail(int64_t tid1, int64_t tid2) const {
    if (thread_nodes.count(llabs(tid1))==0 or thread_nodes.count(llabs(tid2))==0) return {};
    std::vector<int> s;
    std::unordered_set<sgNodeID_t> t2nodes(thread_nodes.at(llabs(tid2)).begin(),thread_nodes.at(llabs(tid2)).end());
    for (auto &nid:thread_nodes.at(llabs(tid1))) {
        s.emplace_back(t2nodes.count(nid));
    }
    return s;
}

void RTGPartition::classify_neighbours(int skip_nodes, bool discard_invalid) {
#pragma omp parallel for
    for (auto b=0;b<thread_neighbours.bucket_count();++b){
        for (auto bi=thread_neighbours.begin(b);bi!=thread_neighbours.end(b);++bi){
            auto &tid=bi->first;
            auto &nns=bi->second;
            std::vector<sgNodeID_t> valid_nodes;
            std::vector<ThreadOverlapType> valid_nodes_types;
            for (auto otid:nns){
                auto oc=rtg.classify_thread_overlap(tid,otid,skip_nodes);
                if (not discard_invalid or oc!=ThreadOverlapType::invalid) {
                    valid_nodes.emplace_back(otid);
                    valid_nodes_types.emplace_back(oc);
                }
            }
#pragma omp critical
            {
                std::swap(nns,valid_nodes);
                thread_neighbours_types[tid]=valid_nodes_types;
            }
        }
    }
}

bool RTGPartition::thread_available(uint64_t tid, int min_nodes, float max_classified_nodes_perc) {
    auto tc=thread_class[tid];
    if (tc!=0) return false;
    auto &tns=thread_nodes[tid];
    if (tns.size()<min_nodes) return false;
    auto skip_countdown= ceil(max_classified_nodes_perc*tns.size());
    for (auto nid:tns) {
        if (node_class[nid]!=0) {
            --skip_countdown;
            if (skip_countdown==0) break;
        }
    }
    return skip_countdown!=0;

}