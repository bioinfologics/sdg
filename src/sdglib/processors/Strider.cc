//
// Created by Bernardo Clavijo (EI) on 01/05/2020.
//

#include <sdglib/views/NodeView.hpp>
#include "Strider.hpp"
inline uint32_t get_votes(std::unordered_map<sgNodeID_t,uint32_t> votes, sgNodeID_t node){
    auto it=votes.find(node);
    if (it==votes.end()) return 0;
    else return it->second;
}

SequenceDistanceGraphPath Strider::walk_out(sgNodeID_t n) {
    SequenceDistanceGraphPath p(ws.sdg);
    std::unordered_map<sgNodeID_t,uint32_t> votes;
    for (auto &p:ws.paired_reads_datastores[0].mapper.all_paths_fw(n)) {
        for (auto &n:p) ++votes[n];
    }

    p.nodes.emplace_back(n);
    while (true){
        auto nv=ws.sdg.get_nodeview(p.nodes.back());
        uint32_t total_votes=0;
        sgNodeID_t winner=0;
        uint32_t winner_votes=0;
        for (const auto & l: nv.next()){
            auto v=get_votes(votes,l.node().node_id());
            total_votes+=v;
            if (v==winner_votes) winner=0;
            if (v>winner_votes) {
                winner=l.node().node_id();
                winner_votes=v;
            }
        }
        if (winner==0 or total_votes<2 or winner_votes<total_votes*.75 or
            std::find(p.nodes.begin(),p.nodes.end(),winner)!=p.nodes.end()) break;
        p.nodes.emplace_back(winner);
    }
    return p;
}

SequenceDistanceGraphPath Strider::walk_out_in_order(sgNodeID_t n,bool use_pair, bool collapse_pair) {
    auto read_paths=ws.paired_reads_datastores[0].mapper.all_paths_fw(n);
    std::vector<int> next_node(read_paths.size());
    std::vector<bool> past_jump(read_paths.size());

    SequenceDistanceGraphPath p(ws.sdg);


    p.nodes.emplace_back(n);
    while (true){
        //each read votes for its next node, votes before the jump count double
        std::unordered_map<sgNodeID_t,uint32_t> votes;
        for (auto i=0;i<read_paths.size();++i) {
            if (next_node[i] != -1) { //i.e. if the read has not finished already :D
                if (read_paths[i][next_node[i]] == 0) {
                    ++next_node[i];
                    past_jump[i] = true;
                }
                if (not past_jump[i]) votes[read_paths[i][next_node[i]]] += 2;
                else ++votes[read_paths[i][next_node[i]]];
            }

        }

        auto nv=ws.sdg.get_nodeview(p.nodes.back());
        uint32_t total_votes=0;
        sgNodeID_t winner=0;
        uint32_t winner_votes=0;
        for (const auto & l: nv.next()){
            auto v=get_votes(votes,l.node().node_id());
            total_votes+=v;
            if (v==winner_votes) winner=0;
            if (v>winner_votes) {
                winner=l.node().node_id();
                winner_votes=v;
            }
        }
        if (winner==0 or total_votes<2 or winner_votes<total_votes*.75) break;
        p.nodes.emplace_back(winner);
        for (auto i=0;i<read_paths.size();++i) {
            if (read_paths[i][next_node[i]]==winner) ++next_node[i];
        }
    }
    return p;
}