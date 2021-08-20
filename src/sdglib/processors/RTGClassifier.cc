//
// Created by Bernardo Clavijo (EI) on 06/08/2021.
//

#include "RTGClassifier.hpp"

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
RTGClassifier::RTGClassifier(const ReadThreadsGraph &rtg, int min_node_threads, float node_min_percentage, int thread_p,
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
        auto t=rtg.get_thread(tc.first);
        thread_nodes[tc.first].reserve(t.size());
        for (auto &np:t) thread_nodes[tc.first].emplace_back(llabs(np.node));
    }
    sdglib::OutputLog()<<"done with "<<node_class.size()<<" nodes in "<<thread_class.size()<<" threads"<<std::endl;
}


int64_t RTGClassifier::get_node_class(sgNodeID_t nid) {
    if (node_class.count(nid)==0) return -1;
    return node_class[nid];
}
int64_t RTGClassifier::compute_node_class(sgNodeID_t nid) {
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
    //TODO: check for a second?
    if (highest_votes>=min_node_threads and highest_votes>=node_min_percentage*threads.size())
        return highest_votes_class;
    return 0;
}

bool RTGClassifier::switch_node_class(sgNodeID_t nid, int64_t c) {
    if (node_class.count(nid)==0 or node_class[nid]==c) return false;
    node_class[nid]=c;
    for (auto &tid:node_threads[nid]) threads_to_evaluate.insert(tid);
    return true;
}

int64_t RTGClassifier::get_thread_class(int64_t tid) {
    if (thread_class.count(tid)==0) return -1;
    return thread_class[tid];
}

int64_t RTGClassifier::compute_thread_class(int64_t tid) {
    std::map<int64_t, int64_t> class_votes;
    int64_t highest_votes=0;
    int64_t highest_votes_class=0;
    int64_t v;
    int64_t c;
    auto &nodes=thread_nodes[tid];
    std::vector<int64_t> nc;
    nc.reserve(nodes.size());
    for (auto nid:nodes) nc.emplace_back(node_class[nid]);
    for (auto c:nc) {
        if (c==0) continue;
        v = ++class_votes[c];
        if (highest_votes < v) {
            highest_votes = v;
            highest_votes_class = c;
        }
    }
    if (2*highest_votes>nodes.size()) //condition 1: 50+% nodes in class
        return highest_votes_class;
    if (highest_votes>=thread_p) { //condition 2: only one class with p in q
        auto p_in_q=find_all_p_in_q(thread_p,thread_q,nc,true);
        if (p_in_q.size()==1) return p_in_q[0];
    }
    return 0;
}

bool RTGClassifier::switch_thread_class(int64_t tid, int64_t c) {
    if (thread_class.count(tid)==0 or thread_class[tid]==c) return false;
    thread_class[tid]=c;
    for (auto &nid:thread_nodes[tid]) nodes_to_evaluate.insert(nid);
    return true;
}

uint64_t RTGClassifier::propagate(uint64_t steps, bool verbose) {
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

uint64_t RTGClassifier::thread_propagate(uint64_t steps, float vote_perc, int max_noise, bool verbose) {

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

std::unordered_map<std::pair<int64_t, int64_t>, std::vector<int64_t>> RTGClassifier::find_class_bridges(int p, int q) {
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

void RTGClassifier::compute_thread_intersections(int min_threads, int max_threads) {
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

void RTGClassifier::compute_thread_neighbours(int min_shared) {
    thread_neighbours.clear();
    for (auto &ti:thread_intersections) {
        if (ti.second>=min_shared){
            thread_neighbours[ti.first.first].emplace_back(ti.first.second);
            thread_neighbours[ti.first.second].emplace_back(ti.first.first);
        }
    }
}

void RTGClassifier::compute_thread_neighbours_p_q(int p, int q) {
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

std::vector<int64_t> RTGClassifier::get_thread_neighbours(int64_t tid) const {
    if (thread_neighbours.count(llabs(tid))==0) return {};
    return thread_neighbours.at(llabs(tid));
}

void RTGClassifier::reset(int _min_node_threads, float _node_min_percentage, int _thread_p, int _thread_q) {
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

int64_t RTGClassifier::get_thread_intersection(int64_t tid1, int64_t tid2) const {
    std::pair<int64_t, int64_t> p;
    tid1 = llabs(tid1);
    tid2 = llabs(tid2);
    if (tid1<tid2) p={tid1,tid2};
    else p={tid2,tid1};
    if (thread_intersections.count(p)) return thread_intersections.at(p);
    return 0;
}

std::vector<int> RTGClassifier::thread_shared_detail(int64_t tid1, int64_t tid2) const {
    if (thread_nodes.count(llabs(tid1))==0 or thread_nodes.count(llabs(tid2))==0) return {};
    std::vector<int> s;
    std::unordered_set<sgNodeID_t> t2nodes(thread_nodes.at(llabs(tid2)).begin(),thread_nodes.at(llabs(tid2)).end());
    for (auto &nid:thread_nodes.at(llabs(tid1))) {
        s.emplace_back(t2nodes.count(nid));
    }
    return s;
}