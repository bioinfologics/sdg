//
// Created by Bernardo Clavijo (EI) on 01/05/2020.
//

#include <sdglib/views/NodeView.hpp>
#include "Strider.hpp"

const std::string Strider::logo="                    80CLfft11G\n                   01;;1111tf8\n                   0;:1CCCCCC0\n                    L:tCCCCCCC8\n                     CifCCCCCC8\n                      8GCCCCC0\n                   8888GLffLG8        088\n                80GGCCCL11111fG8     0LL0\n               0CCCCCCCL11111tLG  80GCG0\n                8GCCC00f111111tLC0GCCG8\n                 8GCC8G1111111tLCCCG8\n                  0CCf1111111CGCG8\n                   0Cf1iii11L\n                   G1CLt:,:;L0\n                   C:;i;,,,,,:iL8\n                   t,,,,,,,,,,ifC0\n      00          0:,,,,,i;:ifCCCC0\n    8L11L0888     G:,,,,i88G8 8CCCG\n  0t1tfCCCCCGGGGGGLf1;:C     8CCCC8\n 0i1C08 800GCCCCCCCCCG0      8CCCC0\n LG        88000GG08         0CCC0\n                              CCC8\n                              GCC8\n                             8fffL08\n                             0t1t1ttfC\n                              8800G08\n";

Strider::Strider(WorkSpace & _ws):ws(_ws){
};

void Strider::dump(std::string filename) {
    std::ofstream ofs(filename);
    sdglib::write_flat_vectorvector(ofs,routes_fw);
    sdglib::write_flat_vectorvector(ofs,routes_bw);
    sdglib::write_flat_vectorvector(ofs,links_fw);
    sdglib::write_flat_vectorvector(ofs,links_bw);
    sdglib::write_bool_vector(ofs,is_anchor);
}

void Strider::load(std::string filename) {
    std::ifstream ifs(filename);
    sdglib::read_flat_vectorvector(ifs,routes_fw);
    sdglib::read_flat_vectorvector(ifs,routes_bw);
    sdglib::read_flat_vectorvector(ifs,links_fw);
    sdglib::read_flat_vectorvector(ifs,links_bw);
    sdglib::read_bool_vector(ifs,is_anchor);
}

inline uint32_t get_votes(std::unordered_map<sgNodeID_t,uint32_t> votes, sgNodeID_t node){
    auto it=votes.find(node);
    if (it==votes.end()) return 0;
    else return it->second;
}

SequenceDistanceGraphPath Strider::stride_out(sgNodeID_t n) {
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

SequenceDistanceGraphPath Strider::stride_out_in_order(sgNodeID_t n,bool use_pair, bool collapse_pair, bool verbose) {
    std::vector<std::vector<sgNodeID_t>> read_paths;
    for (auto prds:paired_datastores) {
        auto new_read_paths = prds->mapper.all_paths_fw(n, use_pair, collapse_pair);
        read_paths.insert(read_paths.end(),new_read_paths.begin(),new_read_paths.end());
    }
    for (auto lrr:long_recruiters) {
        auto new_read_paths = lrr->all_paths_fw(n);
        read_paths.insert(read_paths.end(),new_read_paths.begin(),new_read_paths.end());
    }
    //TODO: improve paths -> add missing small connecting nodes.
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
        if (winner==0 or total_votes<2 or winner_votes<total_votes*.75) {
            if (verbose) {
                std::cout << "  can't find winner, options are:";
                for (const auto & l: nv.next()) std::cout<<" "<<l;
                std::cout<<std::endl;
                std::cout << "  but votes are:";
                for (const auto & v: votes) std::cout<<"    "<<v.first<<": "<<v.second<<std::endl;
            }
            if (!experimental_striding) break;
            //allow to skip a single node (i.e. add votes to the next node to this node's votes too:
            for (const auto & l: nv.next()){
                auto v=get_votes(votes,l.node().node_id());
                for (auto nl: l.node().next()) v+=get_votes(votes,nl.node().node_id());
                total_votes+=v;
                if (v==winner_votes) winner=0;
                if (v>winner_votes) {
                    winner=l.node().node_id();
                    winner_votes=v;
                }
            }
            if (winner==0 or total_votes<2 or winner_votes<total_votes*.75) break;
        }
        p.nodes.emplace_back(winner);
        for (auto i=0;i<read_paths.size();++i) {
            if (next_node[i]!=-1 and read_paths[i][next_node[i]]==winner) {
                if (++next_node[i]==read_paths[i].size()) next_node[i]=-1;
            }
        }
    }
    return p;
}

void Strider::stride_from_anchors(uint32_t min_size, float min_kci, float max_kci) {
    std::cout<<std::endl<<logo<<std::endl;
    sdglib::OutputLog()<<"Gone striding..."<<std::endl;
    routes_fw.clear();
    routes_bw.clear();
    is_anchor.clear();
    routes_fw.resize(ws.sdg.nodes.size());
    routes_bw.resize(ws.sdg.nodes.size());
    is_anchor.resize(ws.sdg.nodes.size());
    if (min_size==0) min_size=1;//a min_size of 0 would get deleted nodes and such
    uint64_t anchors=0,found_fw=0,found_bw=0;
#pragma omp parallel for reduction(+:anchors) reduction(+:found_fw) reduction(+:found_bw)
    for (auto nid=1;nid<ws.sdg.nodes.size();++nid) {
        if (ws.sdg.get_node_size(nid)<min_size) continue;
        auto nv=ws.sdg.get_nodeview(nid);
        if (nv.kci()<min_kci or nv.kci()>max_kci) continue;
        ++anchors;
        is_anchor[nid]=true;
        routes_fw[nid]=stride_out_in_order(nid).nodes;
        if (routes_fw[nid].size()>1) ++found_fw;
        routes_bw[nid]=stride_out_in_order(-nid).nodes;
        if (routes_bw[nid].size()>1) ++found_bw;
    }
    sdglib::OutputLog()<<"Strider found "<<found_fw<<" forward and "<<found_bw<<" backward routes from "<< anchors << " anchors"<<std::endl;
}

std::vector<Link> Strider::link_out_by_lr(sgNodeID_t n,int d, int min_reads, int group_size, int small_node_size, float candidate_percentaje, float first_percentaje, bool verbose) {
    PerfectMatchesMergeSorter pmms(ws);
    for (auto &llr:long_recruiters) {
        pmms.init_from_node(n, *llr, min_reads, group_size, small_node_size);
        std::vector<int32_t> prev_status;
        for (pmms.find_next_node(d, candidate_percentaje, first_percentaje, verbose); pmms.next_node!=0 and pmms.read_next_match!=prev_status; pmms.find_next_node(d, candidate_percentaje, first_percentaje, verbose)){
            prev_status=pmms.read_next_match;//This is a cheap exit to avoid when pmms get stuck in the last node
            pmms.advance_reads_to_node();
            pmms.advance_reads_through_node();
        }
    }
    std::vector<bool> used(pmms.out.size());
    std::vector<Link> links;
    std::vector<int64_t> offsets;
    offsets.reserve(pmms.out.size());
    uint64_t last_nodepos=0;
    for (auto i=0;i<pmms.out.size();++i){
        if (used[i]) continue;
        int skipped=0;
        offsets.clear();
        for (auto j=i;j<pmms.out.size();++j){
            if (pmms.out[j].node==pmms.out[i].node and pmms.out[j].node_position>=last_nodepos){
                offsets.emplace_back(pmms.out[j].read_position-pmms.out[j].node_position);
                last_nodepos=pmms.out[j].node_position;
                used[j]=true;
            }
            else {
                ++skipped;
                if (skipped>=5) break;
            }
        }
        if (offsets.size()>5) {
            //mean for first 10 is used as distance
            //TODO: use mean of last 10 to correct the distance to the next node
            int64_t t=0,c=0;
            for (auto i=0;i<10 and i<offsets.size();++i){
                t+=offsets[i];
                ++c;
            }
            if (links.empty() or links.back().dest != pmms.out[i].node or t/c - links.back().dist > ws.sdg.get_node_size(pmms.out[i].node)*.9)
            links.emplace_back(-n,pmms.out[i].node,t/c);
        }
    }
    return links;
}

void Strider::link_from_anchors(uint32_t min_size, float min_kci, float max_kci, int d, int min_reads, int group_size, int small_node_size, float candidate_percentaje ,float first_percentaje) {
    std::cout<<std::endl<<logo<<std::endl;
    sdglib::OutputLog()<<"Gone linking..."<<std::endl;
    links_fw.clear();
    links_bw.clear();
    is_anchor.clear();
    links_fw.resize(ws.sdg.nodes.size());
    links_bw.resize(ws.sdg.nodes.size());
    is_anchor.resize(ws.sdg.nodes.size());
    if (min_size==0) min_size=1;//a min_size of 0 would get deleted nodes and such
    uint64_t anchors=0,found_fw=0,found_bw=0;
#pragma omp parallel for reduction(+:anchors) reduction(+:found_fw) reduction(+:found_bw)
    for (auto nid=1;nid<ws.sdg.nodes.size();++nid) {
        if (ws.sdg.get_node_size(nid)<min_size) continue;
        auto nv=ws.sdg.get_nodeview(nid);
        if (nv.kci()<min_kci or nv.kci()>max_kci) continue;
        ++anchors;
        is_anchor[nid]=true;
        links_fw[nid]=link_out_by_lr(nid, d, min_reads, group_size, small_node_size, candidate_percentaje, first_percentaje);
        if (!links_fw[nid].empty()) ++found_fw;
        links_bw[nid]=link_out_by_lr(-nid, d, min_reads, group_size, small_node_size, candidate_percentaje, first_percentaje);
        if (!links_bw[nid].empty()) ++found_bw;
    }
    sdglib::OutputLog()<<"Strider found "<<found_fw<<" forward and "<<found_bw<<" backward routes from "<< anchors << " anchors"<<std::endl;
}

void Strider::route_vs_readpaths_stats() {
    uint64_t p_finished80=0;
    uint64_t end_in_tip=0;
    uint64_t reaching_anchor=0;
    uint64_t not_reaching_paths_anchors=0;
    uint64_t anchor_routes=0;
    for (auto fnid=1;fnid<is_anchor.size();++fnid){
        if (not is_anchor[fnid]) continue;
        anchor_routes+=2;
        for (auto nid:{fnid,-fnid}) {
            uint64_t p_finished=0,p_anchor=0;
            std::vector<std::vector<sgNodeID_t>> read_paths;
            //get_readpaths
            for (auto prds:paired_datastores) {
                auto new_read_paths = prds->mapper.all_paths_fw(nid);
                read_paths.insert(read_paths.end(), new_read_paths.begin(), new_read_paths.end());
            }
            for (auto lrr:long_recruiters) {
                auto new_read_paths = lrr->all_paths_fw(nid);
                read_paths.insert(read_paths.end(), new_read_paths.begin(), new_read_paths.end());
            }
            std::vector<sgNodeID_t> route;
            if (nid==fnid) route=routes_fw[fnid];
            else route=routes_bw[fnid];

            //Percentaje of paths with their last node in the route (proxy for "finished reads")
            std::unordered_set<sgNodeID_t> nodes_in_route(route.begin(),route.end());
            for (auto &p:read_paths) if (nodes_in_route.count(p.back())) ++p_finished;
            if (p_finished>=read_paths.size()*.8) ++p_finished80;

            //Percentaje of routes finishing in a tip
            if (ws.sdg.get_nodeview(route.back()).next().size()==0) ++end_in_tip;

            //Percentaje of routes reaching an anchor - reaching an anchor if there's an anchor on 80+% of their paths
            bool anchor_reached=false;
            for (auto &n:route) if (is_anchor[llabs(n)]) {
                if (n==nid) continue;
                ++reaching_anchor;
                anchor_reached=true;
                break;
            }
            if (not anchor_reached) {
                uint64_t anchor_not_in_route = 0;
                for (auto &p:read_paths) {
                    for (auto &n:p) {
                        if (n==nid) continue;
                        if (is_anchor[llabs(n)]) {
                            if (not nodes_in_route.count(n)) ++anchor_not_in_route;
                            break;
                        }
                    }

                }
                if (anchor_not_in_route >= read_paths.size() * .5) ++not_reaching_paths_anchors;
            }
        }
    }

    std::cout<<p_finished80<<" / "<<anchor_routes<<" routes from anchors have 80+% of their read paths's last nodes"<<std::endl;
    std::cout<<end_in_tip<<" routes end in tips"<<std::endl;
    std::cout<<reaching_anchor<<" routes reach other anchors"<<std::endl;
    std::cout<<not_reaching_paths_anchors<<" routes do not reach anchors, but 50+% of their paths do"<<std::endl;
}