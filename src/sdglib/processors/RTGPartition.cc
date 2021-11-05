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

bool RTGPartition::supported_thread(int64_t tid, int min_support) {
    std::map<int64_t, int> thread_nodes_counted;
    bool last_node_unsafe=false;
    auto tnf=rtg.get_thread_nodes(tid);
    for (int i=0;i<tnf.size();++i){
        if (i>tnf.size()/2) {
            //compute "supported threads" in this node
            int good_threads=0;
            for (auto otid:rtg.node_threads(tnf[i],true)) if (thread_nodes_counted[otid]>=min_support) ++good_threads;
            if (good_threads>=min_support){
                last_node_unsafe=false;
            }
            else {
                if (last_node_unsafe) return false;
                last_node_unsafe=true;
            }
        }
        for (auto otid:rtg.node_threads(tnf[i],true)) if (otid!=tid) ++thread_nodes_counted[otid];
    }
    thread_nodes_counted.clear();
    for (int i=tnf.size()-1;i<=0;--i){
        if (i<tnf.size()/2) {
            //compute "supported threads" in this node
            int good_threads=0;
            for (auto otid:rtg.node_threads(tnf[i],true)) if (thread_nodes_counted[otid]>=min_support) ++good_threads;
            if (good_threads>=min_support){
                last_node_unsafe=false;
            }
            else {
                if (last_node_unsafe) return false;
                last_node_unsafe=true;
            }
        }
        for (auto otid:rtg.node_threads(tnf[i],true)) if (otid!=tid) ++thread_nodes_counted[otid];
    }
    return true;
}

void RTGPartition::classify_all_threads(int min_nodes, float max_classified_nodes_perc) {
    auto ut=find_unclassified_threads(min_nodes,max_classified_nodes_perc);
    for (auto i=0;i<ut.size();++i) {
        auto tid=ut[i];
        if (thread_available(tid,min_nodes,max_classified_nodes_perc) and supported_thread(tid)) {
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

void RTGPartition::reset(int _min_node_threads, float _node_min_percentage, int _thread_p, int _thread_q) {
    if (_min_node_threads!=-1) min_node_threads=_min_node_threads;
    if (_node_min_percentage!=-1) node_min_percentage=_node_min_percentage;
    if (_thread_p!=-1) thread_p=_thread_p;
    if (_thread_q!=-1) thread_p=_thread_q;
    for (auto &tc:thread_class) tc.second=0;
    for (auto &nc:node_class) nc.second=0;
    nodes_to_evaluate.clear();
    threads_to_evaluate.clear();
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