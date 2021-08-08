//
// Created by Bernardo Clavijo (EI) on 06/08/2021.
//

#include "RTGClassifier.hpp"

template <typename T>
std::vector<T> find_all_p_in_q(int p, int q, const std::vector<T>& input){
    std::vector<T> v;
    std::unordered_map<T,int64_t> votes;
    for (auto &x:input) ++votes[x];
    std::vector<int64_t> pos;
    for (auto xc:votes) if (xc.second>p) {
        for (auto i=0;i<input.size();++i) {
            auto x=input[i];
            if (x == xc.first) {
                pos.emplace_back(i);
                if (pos.size()>=p and pos.back()-*(pos.end()-p)<=q) {
                    v.emplace_back(x);
                    break;
                }
            }
        }
        pos.clear();
    }
    return v;
};

RTGClassifier::RTGClassifier(const ReadThreadsGraph &rtg, int min_node_threads, float node_min_percentage, int thread_p,
                             int thread_q): rtg(rtg), min_node_threads(min_node_threads),
                             node_min_percentage(node_min_percentage), thread_p(thread_p), thread_q(thread_q) {

    for (auto &nv:rtg.get_all_nodeviews(false,false)) {
        auto nid=nv.node_id();
        auto nts=rtg.node_threads(nid);
        node_threads[nid].reserve(nts.size());
        node_class[nid]=0;
        for (auto &tid : nts) {
            node_threads[nid].emplace_back(tid);
            thread_class[tid]=0;
        }
    }
    for (auto & tc:thread_class){
        auto t=rtg.get_thread(tc.first);
        thread_nodes[tc.first].reserve(t.size());
        for (auto &np:t) thread_nodes[tc.first].emplace_back(llabs(np.node));
    }
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
    if (2*highest_votes_class>threads.size() or (highest_votes>=min_node_threads and highest_votes>=node_min_percentage*threads.size()))
        return highest_votes_class;
    return 0;
}

bool RTGClassifier::switch_node_class(sgNodeID_t nid, int64_t c) {
    if (node_class[nid]==c) return false;
    node_class[nid]=c;
    for (auto &tid:node_threads[nid]) threads_to_evaluate.insert(tid);
    return true;
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
    if (2*highest_votes_class>nodes.size()) //condition 1: 50+% nodes in class
        return highest_votes_class;
    if (highest_votes_class>=thread_p) { //condition 2: only one class with p in q
        auto p_in_q=find_all_p_in_q(thread_p,thread_q,nc);
        if (p_in_q.size()==1) return p_in_q[0];
    }
    return 0;
}

bool RTGClassifier::switch_thread_class(int64_t tid, int64_t c) {
    if (thread_class[tid]==c) return false;
    thread_class[tid]=c;
    for (auto &nid:thread_nodes[tid]) nodes_to_evaluate.insert(nid);
    return true;
}

uint64_t RTGClassifier::propagate(uint64_t steps, bool verbose) {
    uint64_t switched;
    while (steps-- and threads_to_evaluate.size()+nodes_to_evaluate.size()) {
        std::unordered_set<int64_t> otte;
        std::swap(otte,threads_to_evaluate);
        for (auto tid:otte) {
            if(switch_thread_class(tid,compute_thread_class(tid))) ++switched; //switch won't do anything if class doesn't change
        }
        std::unordered_set<sgNodeID_t> onte;
        std::swap(onte,nodes_to_evaluate);
        for (auto nid:onte)
            if (switch_node_class(nid,compute_node_class(nid))) ++switched; //switch won't do anything if class doesn't change
    }
    return switched;
}