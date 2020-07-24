//
// Created by Bernardo Clavijo (EI) on 08/11/2019.
//

#include "GraphContigger.hpp"
#include "GraphEditor.hpp"
#include <sdglib/views/NodeView.hpp>


void GraphContigger::pop_bubbles(const PairedReadsDatastore &prds, int bubble_size, int min_support, int max_noise,
                                    float snr) {
    std::set<sgNodeID_t> to_delete;
    for (auto &nv:ws.sdg.get_all_nodeviews()) {
        if (nv.size() > bubble_size) continue;
        auto p = nv.parallels();
        if (p.size() != 1) continue;
        auto noise = prds.mapper.paths_in_node[nv.node_id()].size();
        if (noise > max_noise) continue;
        auto support = prds.mapper.paths_in_node[llabs(p[0].node_id())].size();
        if (support < min_support or (noise > 0 and (float) support / noise < snr)) continue;
        to_delete.insert(nv.node_id());
    }
    for (auto n:to_delete) ws.sdg.remove_node(n);
    ws.sdg.join_all_unitigs();
}

void GraphContigger::solve_canonical_repeats_with_single_paths(const PairedReadsDatastore & prds,int min_support, int max_noise, float snr, bool join_unitigs, bool dry_run, bool verbose) {
    GraphEditor ge(ws);
    uint64_t repeats=0;
    uint64_t solved_repeats=0;
    for (auto & nv :ws.sdg.get_all_nodeviews()){
        auto nvp=nv.prev();
        auto nvn=nv.next();
        if (nvp.size()!=2 or nvn.size()!=2) continue;//XXX: it should check the nodes in and out have only one connection to this node too!!!!
        if (nvp[0].node().next().size()!=1
            or nvp[1].node().next().size()!=1
            or nvn[0].node().prev().size()!=1
            or nvn[1].node().prev().size()!=1)
            continue;
        //sdglib::OutputLog()<<"repeat on node "<<nv<<std::endl;
        ++repeats;
        std::unordered_map<std::pair<sgNodeID_t,sgNodeID_t>,uint64_t> through;
        for (auto &pid:prds.mapper.paths_in_node[nv.node_id()]){
            //detail_log<<"checking path "<<pid<<std::endl;
            auto p=prds.mapper.read_paths[llabs(pid)].path;
            if (pid<0) {
                auto old_p=p;
                p.clear();
                for (auto on=old_p.rbegin();on!=old_p.rend();++on) p.emplace_back(-(*on));
            }
            if (std::count(p.cbegin(),p.cend(),nv.node_id())!=1) continue;
            auto pnp=std::find(p.cbegin(),p.cend(),nv.node_id());
            if (pnp==p.cbegin() or pnp>=p.end()-1) continue;
            std::pair<sgNodeID_t,sgNodeID_t> t(*(pnp-1),*(pnp+1));
            if (not through.count(t))through[t]=0;
            ++through[t];
        }
        //sdglib::OutputLog()<<"through collection created "<<nv<<std::endl;
        //Check there is 2 winners with adequate support, snr is ok, and no noise has support
        if (through.size()<2) continue;
        std::vector<std::pair<uint64_t,std::pair<sgNodeID_t,sgNodeID_t>>> through_ranking;
        for (auto &t:through) through_ranking.emplace_back(t.second,t.first);
        std::sort(through_ranking.rbegin(),through_ranking.rend());

        //sdglib::OutputLog()<<"through ranking created "<<nv<<std::endl;
        //for (auto &tr:through_ranking) sdglib::OutputLog()<<" "<<tr.second.first<<" <-> "<<tr.second.second<<"  x"<<tr.first<<std::endl;
        if (through_ranking[1].first<min_support) continue;
        //sdglib::OutputLog()<<"."<<std::endl;
        if (through_ranking.size()>2 and (through_ranking[2].first>max_noise or through_ranking[2].first*snr>=through_ranking[1].first)) continue;
        //sdglib::OutputLog()<<"."<<std::endl;
        if (through_ranking[0].second.first==through_ranking[1].second.first or through_ranking[0].second.second==through_ranking[1].second.second) continue;
        //sdglib::OutputLog()<<"."<<std::endl;
        //Check the winners are the in prev/next
        if ( not ( (through_ranking[0].second.first==nvp[0].node().node_id() and through_ranking[1].second.first==nvp[1].node().node_id()) or
                (through_ranking[1].second.first==nvp[0].node().node_id() and through_ranking[0].second.first==nvp[1].node().node_id()))) continue;
        //sdglib::OutputLog()<<"."<<std::endl;
        if ( not ( (through_ranking[0].second.second==nvn[0].node().node_id() and through_ranking[1].second.second==nvn[1].node().node_id()) or
                   (through_ranking[1].second.second==nvn[0].node().node_id() and through_ranking[0].second.second==nvn[1].node().node_id()))) continue;
        //add expansion to ge

        //sdglib::OutputLog()<<"queueing expansion "<<nv<<std::endl;
        ge.queue_node_expansion(nv.node_id(),{{-through_ranking[0].second.first,through_ranking[0].second.second},{-through_ranking[1].second.first,through_ranking[1].second.second}});
        ++solved_repeats;
        if (verbose){
            std::cout << "---" << std::endl;
            std::cout << "Solving: " << nv.node_id() << std::endl;
            std::cout << "seq"<<llabs(through_ranking[0].second.first) << ",seq"<<llabs(nv.node_id())<<",seq"<<llabs(through_ranking[0].second.second) << std::endl;
            std::cout << "seq"<<llabs(through_ranking[1].second.first) << ",seq"<<llabs(nv.node_id())<<",seq"<<llabs(through_ranking[1].second.second) << std::endl;
        }
    }
    sdglib::OutputLog()<<"Contigger repeat_expansion: "<<solved_repeats<<" / "<<repeats<<std::endl;
    if (!dry_run){
        ge.apply_all();
    }

    if (join_unitigs){
        ws.sdg.join_all_unitigs();
    }

}

void GraphContigger::solve_canonical_repeats_with_paired_paths(const PairedReadsDatastore & prds,int min_support, int max_noise, float snr, bool join_unitigs, bool dry_run=false, bool verbose=false) {
    GraphEditor ge(ws);
    uint64_t repeats=0;
    uint64_t solved_repeats=0;
    for (auto & nv :ws.sdg.get_all_nodeviews()) {
        auto nvp = nv.prev();
        auto nvn = nv.next();
        std::vector<std::pair<sgNodeID_t, sgNodeID_t >> through;
        bool non_canonical = false;
        if (nvp.size()>1 and nvp.size() == nvn.size()) {
            bool stop=false;
            for (auto &l:nvp) if (l.node().next().size()>1) stop=true;
            for (auto &l:nvn) if (l.node().prev().size()>1) stop=true;
            if (stop) continue;
            std::unordered_map<std::pair<sgNodeID_t, sgNodeID_t>, uint64_t> through_counts;
            for (auto &p:nvp) {
                auto pn = p.node().node_id();
                for (auto &pid:prds.mapper.paths_in_node[llabs(pn)]) {
                    for (auto &n:nvn) {
                        auto nn = n.node().node_id();
                        if (pn == nn) {
                            non_canonical = true;

                            break;
                        }
                        for (auto &pidn:prds.mapper.paths_in_node[llabs(nn)]) {
                            if ((llabs(pid) + 1) / 2 == (llabs(pidn) + 1) / 2) {
                                through_counts[{pn, nn}] += 1;
                            }
                        }
                    }
                    if (non_canonical) break;
                }
                if (non_canonical) break;
            }
            //sdglib::OutputLog() << "repeat on node " << nv << std::endl;
            //for (auto &tc:through_counts) {
            //    sdglib::OutputLog() << " " << tc.first.first << "->" << tc.first.second << ":" << tc.second<<std::endl;
            //}


            std::unordered_map<sgNodeID_t, sgNodeID_t> in_best_through;
            std::unordered_map<sgNodeID_t, uint64_t > in_best_score;
            std::unordered_map<sgNodeID_t, uint64_t > in_second_score;
            std::unordered_map<sgNodeID_t, sgNodeID_t> out_best_through;
            std::unordered_map<sgNodeID_t, uint64_t > out_best_score;
            std::unordered_map<sgNodeID_t, uint64_t > out_second_score;
            for (auto &tc:through_counts) {
                if (tc.second>in_best_score[tc.first.first]){
                    in_second_score[tc.first.first]=in_best_score[tc.first.first];
                    in_best_score[tc.first.first]=tc.second;
                    in_best_through[tc.first.first]=tc.first.second;
                } else if (tc.second>in_second_score[tc.first.first]) in_second_score[tc.first.first]=tc.second;
                if (tc.second>out_best_score[tc.first.second]){
                    out_second_score[tc.first.second]=out_best_score[tc.first.second];
                    out_best_score[tc.first.second]=tc.second;
                    out_best_through[tc.first.second]=tc.first.second;
                } else if (tc.second>out_second_score[tc.first.second]) out_second_score[tc.first.second]=tc.second;
            }
            //Condition 1: min_support -> for every in and out check the best support is >= min
            //Conditions 2&3: snr -> for every in and out check the ratio of the second-best vs. the best is >= snr and support<=max_noise
            if (in_best_score.size()!=nvp.size() or out_best_score.size()!=nvn.size()) continue;
            bool cond_fail=false;
            for (auto &best:in_best_score){
                auto &second=in_second_score[best.first];
                if (second>max_noise or second*snr>best.second or best.first<min_support) cond_fail=true;
            }
            for (auto &best:out_best_score){
                auto &second=out_second_score[best.first];
                if (second>max_noise or second*snr>best.second or best.first<min_support) cond_fail=true;
            }
            if (cond_fail) continue;
            std::vector<std::pair<sgNodeID_t ,sgNodeID_t >> paths_through;
            for (auto &best:in_best_through) paths_through.emplace_back(-best.first,best.second);
            ge.queue_node_expansion(nv.node_id(),paths_through);
            ++solved_repeats;

            if (verbose){
                std::cout << "---" << std::endl;
                std::cout << "Solving: " << nv.node_id() << std::endl;
                std::cout << "seq"<<paths_through[0].first << ",seq"<<llabs(nv.node_id())<<",seq"<<paths_through[0].second << std::endl;
                std::cout << "seq"<<paths_through[1].first << ",seq"<<llabs(nv.node_id())<<",seq"<<paths_through[1].second << std::endl;
            }
            //C
        }
    }
    sdglib::OutputLog()<<"Contigger repeat_expansion: "<<solved_repeats<<" / "<<repeats<<std::endl;
    if (!dry_run){
        ge.apply_all();
    }

    if (join_unitigs){
        ws.sdg.join_all_unitigs();
    }
}

void GraphContigger::solve_canonical_repeats_with_long_reads(const LongReadsRecruiter &lrr, float max_side_kci, int min_support, int max_noise, float snr) {
    GraphEditor ge(ws);
    uint64_t repeats=0;
    uint64_t solved_repeats=0;
    for (auto & nv :ws.sdg.get_all_nodeviews()){
        auto nvp=nv.prev();
        auto nvn=nv.next();
        if (nvp.size()!=2 or nvn.size()!=2) continue;//XXX: it should check the nodes in and out have only one connection to this node too!!!!
        if (nvp[0].node().next().size()!=1
            or nvp[1].node().next().size()!=1
            or nvn[0].node().prev().size()!=1
            or nvn[1].node().prev().size()!=1)
            continue;
        //sdglib::OutputLog()<<"repeat on node "<<nv<<std::endl;

        sgNodeID_t r=nv.node_id(),pA=nvp[0].node().node_id(), pB=nvp[1].node().node_id(), nA=nvn[0].node().node_id(), nB=nvn[1].node().node_id();
        auto pAkci=ws.sdg.get_nodeview(pA).kci();
        if (pAkci==-1 or pAkci>max_side_kci) continue;
        auto pBkci=ws.sdg.get_nodeview(pB).kci();
        if (pBkci==-1 or pBkci>max_side_kci) continue;
        auto nAkci=ws.sdg.get_nodeview(nA).kci();
        if (nAkci==-1 or nAkci>max_side_kci) continue;
        auto nBkci=ws.sdg.get_nodeview(nB).kci();
        if (nBkci==-1 or nBkci>max_side_kci) continue;
        //sdglib::OutputLog()<<"repeat on node "<<nv<<" passed kci conditions, voting"<<std::endl;
        ++repeats;
        //Count threads doing pA->nA, pA->nB and nA->other (TODO: make the conditions more strict!)
        uint64_t vAA=0,vAB=0,vAo=0;
        //sdglib::OutputLog()<<"voting from pA="<<pA<<std::endl;
        if (r==115695) std::cout<<"r="<<r<<" pA="<<pA<<" pB="<<pB<<" nA="<<nA<<" nB="<<nB<<std::endl;
        for (auto tid:lrr.node_threads[llabs(pA)]){
            auto &t=lrr.read_threads[tid];
            if (r==115695) {
                std::cout<<"thread["<<tid<<"] : ";
                for (auto &tm:t) std::cout<<" "<<tm.node;
                std::cout<<std::endl;
            }
            for (auto i=0;i<t.size();++i){
                if (t[i].node==pA and i<t.size()-1){
                    //std::cout<<"found "<<pA<<std::endl;
                    if (t[i+1].node==nA or (i<t.size()-2 and t[i+1].node==r and t[i+2].node==nA)) ++vAA;
                    else if (t[i+1].node==nB or (i<t.size()-2 and t[i+1].node==r and t[i+2].node==nB)) ++vAB;
                    else ++vAo;
                    break;
                }
                if (t[i].node==-pA and i>0){
                    //std::cout<<"found "<<-pA<<std::endl;
                    if (t[i-1].node==-nA or (i>1 and t[i-1].node==-r and t[i-2].node==-nA)) ++vAA;
                    else if (t[i-1].node==-nB or (i>1 and t[i-1].node==-r and t[i-2].node==-nB)) ++vAB;
                    else ++vAo;
                    break;
                }
            }
            if (r==115695) std::cout<<"vAA="<<vAA<<" vAB="<<vAB<<" vAo="<<vAo<<std::endl;
        }
        //Count threads doing B->A, B->B and B->other
        uint64_t vBA=0,vBB=0,vBo=0;
        //sdglib::OutputLog()<<"voting from pB="<<pB<<std::endl;
        for (auto tid:lrr.node_threads[llabs(pB)]){
            auto &t=lrr.read_threads[tid];
            if (r==115695) {
                std::cout<<"thread["<<tid<<"] : ";
                for (auto &tm:t) std::cout<<" "<<tm.node;
                std::cout<<std::endl;
            }
            for (auto i=0;i<t.size();++i){
                if (t[i].node==pB and i<t.size()-1){
                    //std::cout<<"found "<<pB<<std::endl;
                    if (t[i+1].node==nA or (i<t.size()-2 and t[i+1].node==r and t[i+2].node==nA)) ++vBA;
                    else if (t[i+1].node==nB or (i<t.size()-2 and t[i+1].node==r and t[i+2].node==nB)) ++vBB;
                    //else ++vBo;
                    break;
                }
                if (t[i].node==-pB and i>0){
                    //std::cout<<"found "<<-pB<<std::endl;
                    if (t[i-1].node==-nA or (i>1 and t[i-1].node==-r and t[i-2].node==-nA)) ++vBA;
                    else if (t[i-1].node==-nB or (i>1 and t[i-1].node==-r and t[i-2].node==-nB)) ++vBB;
                    //else ++vBo;
                    break;
                }
            }
            if (r==115695) std::cout<<"vBA="<<vBA<<" vBB="<<vBB<<" vBo="<<vBo<<std::endl;
        }
        //sdglib::OutputLog()<<"repeat on node "<<nv<<" voted, evaluating: vAA="<<vAA<<" vAB="<<vAB<<" vAo="<<vAo<<" vBA="<<vBA<<" vBB="<<vBB<<" vBo="<<vBo<<std::endl;
        if (vAA>min_support and vBB>min_support and vAA/snr>=vAB+vAo and vBB/snr>=vBA+vBo) {
            ge.queue_node_expansion(r, {{-pA,nA},{-pB,nB}});
            ++solved_repeats;
        }
        if (vAB>min_support and vBA>min_support and vAB/snr>=vAA+vAo and vBA/snr>=vBB+vBo) {
            ge.queue_node_expansion(r, {{-pA,nB},{-pB,nA}});
            ++solved_repeats;
        }
    }
    sdglib::OutputLog()<<"Contigger LRR repeat_expansion: "<<solved_repeats<<" / "<<repeats<<std::endl;
    ge.apply_all();
    ws.sdg.join_all_unitigs();
}

void GraphContigger::reconnect_tips(const PairedReadsDatastore & prds, int min_support) {
    //first mark all nodes that are tips
    std::vector<bool> tip_fw(ws.sdg.nodes.size());
    std::vector<bool> tip_bw(ws.sdg.nodes.size());

    for (auto &nv:ws.sdg.get_all_nodeviews()){
        if (nv.next().empty()) tip_fw[nv.node_id()]=true;
        if (nv.prev().empty()) tip_bw[nv.node_id()]=true;
    }

    std::unordered_map<std::pair<sgNodeID_t,sgNodeID_t>,uint64_t> pathcount;
    sgNodeID_t t1=0,t2=0;
    bool failed=false;
    for (auto i1=1;i1<prds.mapper.read_paths.size();i1+=2){
        failed=false;
        t1=0;
        t2=0;
        for (auto pit=prds.mapper.read_paths[i1].path.cbegin();pit<prds.mapper.read_paths[i1].path.cend();++pit) {
            if ( (*pit>0 and tip_fw[*pit]) or (*pit<0 and tip_fw[-*pit]) ) {
                if (t1==0 or t1==*pit) t1=*pit;
                else failed=true;
            }
            if ( (*pit>0 and tip_bw[*pit]) or (*pit<0 and tip_bw[-*pit]) ) {
                if (t1!=0 and (t2==0 or t2==*pit)) t2=*pit;
                else failed=true;
            }
        }
        for (auto pit=prds.mapper.read_paths[i1+1].path.crbegin();pit<prds.mapper.read_paths[i1+1].path.crend();++pit) {
            sgNodeID_t n=-*pit;
            if ( (n>0 and tip_fw[n]) or (n<0 and tip_fw[-n]) ) {
                if (t1==0 or t1==n) t1=n;
                else failed=true;
            }
            if ( (n>0 and tip_bw[n]) or (n<0 and tip_bw[-n]) ) {
                if (t1!=0 and (t2==0 or t2==n)) t2=n;
                else failed=true;
            }
        }
        if ( t1!=0 and t2!=0 and not failed ) ++pathcount[{std::min(t1,t2),std::max(t1,t2)}];
    }
    uint64_t r=0;
    std::ofstream tcdf("tip_conn_detail.txt");
    for (auto &pc:pathcount){
        if (pc.second>=min_support) r++;
        tcdf<<pc.first.first<<" "<<pc.first.second<<" "<<pc.second<<std::endl;
    }
    std::cout<<"tip reconnection: "<<r<<" tip pairs connected by at least "<<min_support<<" links"<<std::endl;
}

void GraphContigger::clip_tips(int tip_size, int rounds) {
    for (auto r=0;r<rounds;++r) {
        std::set<sgNodeID_t> to_delete;
        for (auto &nv:ws.sdg.get_all_nodeviews()){
            if (nv.size()<=tip_size and (nv.prev().size()==0 or nv.next().size()==0))
                to_delete.insert(nv.node_id());
        }
        std::cout << "Nodes to delete: " << to_delete.size() << std::endl;
        for (auto n:to_delete) ws.sdg.remove_node(n);
        auto utc = ws.sdg.join_all_unitigs();
        if (to_delete.empty() and utc==0) break;
    }
}

void GraphContigger::remove_small_unconnected(int min_size) {
    for (sgNodeID_t n = 1; n < ws.sdg.nodes.size(); ++n) {
        if (ws.sdg.nodes[n].status == NodeStatus::Deleted) continue;
        if (ws.sdg.nodes[n].sequence.size() >= min_size) continue;
        if (ws.sdg.get_fw_links(n).size()==0 and ws.sdg.get_bw_links(n).size()==0) ws.sdg.remove_node(n);
    }
}

void GraphContigger::extend_to_repeats(int max_size) {
    GraphEditor ge(ws);
    for (auto &fnv:ws.sdg.get_all_nodeviews()){
        for (auto nv:{fnv,fnv.rc()} ){
            if (nv.size()>max_size) continue;
            if (nv.prev().size()<2 or nv.next().size()!=1) continue;
            auto fwn=nv.next()[0].node().node_id();
            bool b=false;
            std::vector<std::pair<sgNodeID_t ,sgNodeID_t >> through;
            for (auto &lp:nv.prev()){
                if (lp.node().next().size()>1 or lp.node().parallels().size()>1){
                    b=true;
                    break;
                }
                through.emplace_back(-lp.node().node_id(),fwn);
            }
            if (b) continue;

            ge.queue_node_expansion(nv.node_id(),through);
        }
    }
    ge.apply_all();
//    def collapse_into_inputs(g,e, max_collapsed_size=500):
//    """Finds all nodes representing only shared sequence at the end of previous nodes, collapses back, only queues operations"""
//    for nvf in g.get_all_nodeviews():
//        for nv in [nvf,nvf.rc()]:
//            if nv.size()>max_collapsed_size: continue #don't expand large nodes
//            if len(nv.prev())<=1 or len(nv.next())!=1: continue
//            b=False
//            for lp in nv.prev():
//                if len(lp.node().next())>1 or len(lp.node().parallels())>1: #a single conflicted prev node skips the expansion, could be solved partially?
//                    b=True
//                    break#don't expand after a bubble
//            if b: continue
//            #TODO: tip check? is it even needed now?
//            #print([[-x.node().node_id(), nv.next()[0].node().node_id()] for x in nv.prev()])
//            e.queue_node_expansion(nv.node_id(),[[-x.node().node_id(), nv.next()[0].node().node_id()] for x in nv.prev()])
}