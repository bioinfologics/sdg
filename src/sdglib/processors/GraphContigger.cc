//
// Created by Bernardo Clavijo (EI) on 08/11/2019.
//

#include "GraphContigger.hpp"
#include "GraphEditor.hpp"
#include <sdglib/views/NodeView.hpp>
#include <sdglib/views/TangleView.hpp>


void GraphContigger::pop_bubbles(const PairedReadsDatastore &prds, int bubble_size, int min_support, int max_noise,
                                    float snr) {

    if (prds.mapper.paths_in_node.size() == 0){
        throw std::runtime_error("Path reads first");
    }
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

bool GraphContigger::solve_bubble(TangleView &t, Strider &s, GraphEditor &ge) {
    // Return false if the bubble is not solved and tue is if't solved and adds an operation to the ge and return the ge response

    std::vector<sgNodeID_t> fapath_all=s.stride_out_in_order(t.frontiers[0].rc().node_id()).nodes;
    int offset = fapath_all.size()<3 ? fapath_all.size() : 3;
    std::vector<sgNodeID_t> fapath (fapath_all.begin(), fapath_all.begin()+offset);

    std::vector<sgNodeID_t> fbpath_all=s.stride_out_in_order(t.frontiers[1].rc().node_id()).nodes;
    offset = fbpath_all.size()<3 ? fbpath_all.size() : 3;
    std::vector<sgNodeID_t> fbpath (fbpath_all.begin(), fbpath_all.begin()+offset);
    std::reverse(fbpath.begin(), fbpath.end());
    std::transform(fbpath.begin(), fbpath.end(), fbpath.begin(), [](sgNodeID_t i){return i*-1;});

    if (fapath.size()<2) fapath=fbpath;
    if (fbpath.size()<2) fbpath=fapath;
    if (fapath.size()<2) return false;

    if (fapath==fbpath and std::abs(t.internals[0].node_id()) == std::abs(fapath[1])){
        return ge.queue_node_deletion(t.internals[1].node_id());
    }
    if (fapath==fbpath and std::abs(t.internals[1].node_id()) == std::abs(fapath[1])){
        return ge.queue_node_deletion(t.internals[0].node_id());
    }
    return false;
}

bool GraphContigger::solve_repeat(TangleView &t, Strider &s, GraphEditor &ge) {
    sgNodeID_t rnid=t.internals[0].node_id();
    std::vector<NodeView> ins;
    std::vector<sgNodeID_t> inids;
    for (const auto& pn: t.internals[0].prev()){
        ins.push_back(pn.node());
        inids.push_back(pn.node().node_id());
    }
    std::vector<NodeView> outs;
    std::vector<sgNodeID_t>outids;
    for (const auto& nn: t.internals[0].next()){
        outs.push_back(nn.node());
        outids.push_back(nn.node().node_id());
    }
    if (ins.size() != outs.size()) return false;

    std::vector<std::vector<sgNodeID_t >> solutions;
    for (const auto& inv: ins){
        auto fpath_all=s.stride_out_in_order(inv.node_id()).nodes;
        int offset = fpath_all.size()<3 ? fpath_all.size() : 3;
        std::vector<sgNodeID_t> fpath(fpath_all.begin(), fpath_all.begin()+offset);
        if (fpath.size()==3 and fpath[1]==rnid and std::find(outids.begin(), outids.end(), fpath[2]) != outids.end()) {
            std::vector<sgNodeID_t> bpath_all=s.stride_out_in_order(-fpath[2]).nodes;
            int offset = bpath_all.size()<3 ? bpath_all.size() : 3;
            std::vector<sgNodeID_t> bpath (bpath_all.begin(), bpath_all.begin()+offset);
            std::reverse(bpath.begin(), bpath.end());
            std::transform(bpath.begin(), bpath.end(), bpath.begin(), [](sgNodeID_t i){return i*-1;});

            if (bpath == fpath) solutions.push_back(fpath);
        }
    }
    if (solutions.size() != ins.size()) return false;
    std::vector<std::pair<sgNodeID_t,sgNodeID_t>> transp_solutions;
    for (const auto& s: solutions){
        transp_solutions.push_back(std::make_pair(-s[0], s[2]));
    }
    ge.queue_node_expansion(rnid,transp_solutions);
    return true;
}

bool GraphContigger::solve_tip(TangleView &t, Strider &s, GraphEditor &ge){
    //just check that each frontier connects to the other one through reads
    auto f0s = s.stride_out_in_order(-t.frontiers[0].node_id()).nodes;
    auto f1s = s.stride_out_in_order(-t.frontiers[1].node_id()).nodes;
    if (f0s.size()>1 and f0s[1]==t.frontiers[1].node_id() and f1s.size()>1 and f1s[1]==t.frontiers[0].node_id()){
        return ge.queue_node_deletion(t.internals[0].node_id());
    }
    return false;
}

bool GraphContigger::solve_unclassified(TangleView &t, Strider &s, GraphEditor &ge){
    if (t.frontiers.size()+t.internals.size()>200) return false;
    std::vector<sgNodeID_t> fnids;
    for (const auto& f: t.frontiers) fnids.push_back(f.node_id());

//    std::unordered_map<sgNodeID_t , std::vector<sgNodeID_t>> sols;
//    for (const auto& fnid: fnids) sols.insert(std::make_pair(fnid, std::vector<sgNodeID_t>()));
    std::unordered_map<sgNodeID_t , std::vector<sgNodeID_t>> fs;
    for (const auto&f: t.frontiers) {
        fs.insert(std::make_pair(-f.node_id(), s.stride_out_in_order(-f.node_id()).nodes));
    }

    std::unordered_map<sgNodeID_t , std::vector<sgNodeID_t>> ns;
    for (const auto&nv: t.internals){

        std::vector<sgNodeID_t> tpath;
        auto stride_back = s.stride_out_in_order(-nv.node_id()).nodes;
        for (auto i=stride_back.rbegin(); i<stride_back.rend()-1; ++i){
            tpath.push_back(*i*-1);
        }
        for (const auto &i: s.stride_out_in_order(nv.node_id()).nodes){
            tpath.push_back(i-1);
        }
        ns.insert(std::make_pair(nv.node_id(), tpath));
    }

    std::vector<std::vector<sgNodeID_t >> sols;
    std::vector<std::vector<sgNodeID_t>> map_values;
    for (const auto &v: fs){map_values.push_back(v.second);}
    for (const auto &v: ns){map_values.push_back(v.second);}


    for (const auto& x: map_values){
        auto sol=end_to_end_solution(x, fnids);
        if (sol.size()>0 and std::find(sols.begin(), sols.end(), sol)==sols.end()){
            sols.push_back(sol);
        }
    }

    std::unordered_set<sgNodeID_t> temp_set;
    for (const auto &x: sols){
        temp_set.emplace(x[0]);
    }
    for (const auto &x: sols){
        temp_set.emplace(-x.back());
    }

    if (sols.size() == fnids.size()/2 and temp_set.size() == fnids.size()){
        for (const auto &x: sols){
            try {
                auto ps=SequenceDistanceGraphPath(ws.sdg,x).sequence();
            } catch (...) {
                return false;
            }
        }
        for (const auto &x: sols) ge.queue_path_detachment(x, true);
        for (const auto &x: t.internals) ge.queue_node_deletion(x.node_id());
        return true;
    }
    return false;
}

std::vector<sgNodeID_t> GraphContigger::end_to_end_solution(std::vector<sgNodeID_t> p, std::vector<sgNodeID_t> fnids){
    for (auto i=0; i<p.size(); ++i){
        if (std::find(fnids.begin(), fnids.end(), -p[i]) != fnids.end()) {
            for (auto j=i+1; j<p.size(); ++j) {
                if (std::find(fnids.begin(), fnids.end(), p[j])!=fnids.end()){
                    std::vector<sgNodeID_t> sol;
                    if (std::abs(p[i])<std::abs(p[j])){
                        for (auto x=i; x<j+1; ++x) sol.push_back(p[x]);
                    } else {
                        for (auto x=i; x<j+1; ++x) sol.push_back(-p[x]);
                        std::reverse(sol.begin(), sol.end());
                    }
                    std::vector<sgNodeID_t> csol;
                    csol.push_back(sol[0]);
                    for (auto x=1; i<sol.size(); ++x){
                        auto n=sol[x];
                        auto fw_reached_nodes_1 = ws.sdg.fw_reached_nodes(csol.back(),1);
                        auto fw_reached_nodes_2 = ws.sdg.fw_reached_nodes(csol.back(),2);
                        if (std::find(fw_reached_nodes_1.begin(), fw_reached_nodes_1.end(), n)!=fw_reached_nodes_1.end()) {
                            csol.push_back(n);
                        } else if (std::find(fw_reached_nodes_2.begin(), fw_reached_nodes_2.end(), n)!=fw_reached_nodes_2.end()) {
                            std::vector<sgNodeID_t> a;
                            for (const auto& skipped: ws.sdg.fw_reached_nodes(csol.back(),1)) {
                                auto fw_skipped = ws.sdg.fw_reached_nodes(skipped,1);
                                if (std::find(fw_skipped.begin(), fw_skipped.end(), n)!=fw_skipped.end()) a.push_back(skipped);
                                if (a.size()==1) {
                                    csol.push_back(a[0]);
                                    csol.push_back(n);
                                } else {
                                    return sol;
                                }
                            }
                        } else {
//                            for (auto &s: sol) {
//                                std::cout << s << ",";
//                            }
//                            std::cout << std::endl;
                            return sol;
                        }
                    }
                    return csol;
                }
            }
            return std::vector<sgNodeID_t>();
        }
    }
    return std::vector<sgNodeID_t>();
}

bool GraphContigger::solve_canonical_repeat(GraphEditor& ge, NodeView &nv, PairedReadsDatastore& peds, int min_support=5, int max_noise=10, int snr=10, bool verbose=false){
    if (verbose) std::cout<< "Solving " << nv.node_id() << ", "<< nv.size() << "bp" << std::endl;
    std::vector<sgNodeID_t > out_nids;
    for (const auto &on: nv.next()) out_nids.push_back(on.node().node_id());

    std::vector<std::pair<sgNodeID_t, sgNodeID_t>> sols;
    for (const auto& lnv: nv.prev()){
        auto pnv = lnv.node();
        std::map<sgNodeID_t,int> c;
        for (const auto& p: peds.mapper.all_paths_fw(pnv.node_id(), false)) {
            if (p[0] != nv.node_id()) {
                c[p[0]]++;
            } else if (p.size() > 1) {
                c[p[1]]++;
            }
        }

        if (verbose) {
            std::cout << "Cuentas de paths" << std::endl;
            for (const auto &cuentas: c) std::cout << cuentas.first <<','<<cuentas.second << std::endl;
        }

        auto t = 0;
        for (const auto &v: c) t+=v.second;
        if (t<min_support){
            if (verbose) std::cout << "FW too few paths! -> " << t << std::endl;
            return false;
        }
        // Pick winner
//        auto winner = std::max_element(c.begin(), c.end(), [](const p1, const p2){return (p1.second<p2.second)});
        int winner=0;
        int votes=0;
        if (verbose) std::cout<<"Counter votes FW -> "<< nv.node_id()<<std::endl;
        for (const auto& counts: c) {
            if (verbose) std::cout<<counts.first<<":"<<counts.second<<std::endl;
            if (votes<counts.second) {
                votes = counts.second;
                winner = counts.first;
            }
        }
        std::cout << "FW Winner: " << winner << ", votes: " << votes << ", total: "<< t << std::endl;
        if (votes<min_support){
            if (verbose) std::cout << "FW winner poorly supported! -> " << winner << ":" << votes << std::endl;
            return false;
        }
        int noise=t-votes;
        if (std::find(out_nids.begin(), out_nids.end(), winner)!=out_nids.end() and noise*snr<=votes and noise<max_noise) {
            sols.push_back(std::make_pair(-pnv.node_id(), winner));
        } else {
            if (verbose) std::cout << "FW Not a clear winner!" << std::endl;
            return false;
        }
    }
    if (verbose) {
        std::cout << "Printing solutions" << std::endl;
        for (const auto &s: sols) {
            std::cout << s.first<<","<<s.second<<std::endl;
        }
    }
    std::unordered_set<sgNodeID_t> unique_ends;
    for (const auto &s: sols) {
        unique_ends.emplace(s.second);
    }
    if (sols.size()<nv.prev().size() or sols.size()>unique_ends.size()) {
        if (verbose) std::cout << "Solutions not match" << std::endl;
        return false;
    }
    std::vector<sgNodeID_t> in_nids;
    for (const auto& x: nv.prev()){
        in_nids.push_back(-x.node().node_id());
    }
    std::vector<std::pair<sgNodeID_t, sgNodeID_t>> sols2;
    for (const auto& lnv: nv.next()){
        auto nnv = lnv.node();
        std::map<sgNodeID_t,int> c;
        for (const auto& p: peds.mapper.all_paths_fw(-nnv.node_id(), false)) {
            if (p[0] != -nv.node_id()) {
                c[p[0]]++;
            } else if (p.size() > 1) {
                c[p[1]]++;
            }
        }
        auto t = 0;
        for (const auto &v: c) t+=v.second;
        if (t<min_support){
            if (verbose) std::cout << "BW too few paths! -> " << t << std::endl;
            return false;
        }
        // Pick winner
//        auto winner = std::max_element(c.begin(), c.end(), [](p1, p2)[return p1.second<p2.second])
        int winner=0;
        int votes=0;
        if (verbose) std::cout<<"Counter votes BW -> "<< nv.node_id()<<std::endl;
        for (const auto& counts: c) {
            if (verbose) std::cout<<counts.first<<":"<<counts.second<<std::endl;
            if (votes<counts.second) {
                votes = counts.second;
                winner = counts.first;
            }
        }
        std::cout << "BW Winner: " << winner << ", votes: " << votes << ", total: "<< t << std::endl;
        if (votes<min_support){
            if (verbose) std::cout << "BW winner poorly supported! ->"<< winner << ":" << votes << std::endl;
            return false;
        }
        int noise=t-votes;
        if (verbose){
            std::cout << "printing in_nids" << std::endl;
            for (const auto &i: in_nids) std::cout << i << std::endl;
        }
        if (std::find(in_nids.begin(), in_nids.end(), winner)!=in_nids.end() and noise*snr<=votes and noise<max_noise) {
            sols2.push_back(std::make_pair(winner, nnv.node_id()));
        } else {
            if (verbose) std::cout << "BW Not a clear winner! -> " << winner << ":" << votes << std::endl;
            return false;
        }
    }
    if (verbose) {
        std::cout << "Printing solutions2" << std::endl;
        for (const auto &s: sols2) {
            std::cout << s.first<<","<<s.second<<std::endl;
        }
    }
    std::sort(sols.begin(), sols.end());
    std::sort(sols2.begin(), sols2.end());
    if (sols!=sols2) {
        if (verbose){
            std::cout << "FW and BW solutions not the same" << std::endl;
            for (auto i=0; i<sols.size();++i){
                std::cout << sols[i].first<<sols[i].second<<sols2[i].first<<sols2[i].second<<std::endl;
            }
        }
        return false;
    }
    ge.queue_node_expansion(nv.node_id(), sols);
    return true;
}

bool GraphContigger::clip_tip(GraphEditor& ge, const NodeView &_nv, PairedReadsDatastore& peds, int min_support, int max_noise, int snr, bool verbose){
    // TODO: small difference vs python version, couldn't find the reason yet
    auto nv (_nv);
    if (nv.prev().size()==0) nv=nv.rc();
    if (nv.prev().size()!=1 or nv.next().size()!=0) return false;
    int ayes=0;
    int nays=0;
    for (const auto &p: peds.mapper.all_paths_fw(nv.prev()[0].node().node_id(), false)){
        if (p[0] == nv.node_id()) {
            ayes++;
        } else {
            nays++;
        }
    }
    if (ayes<=max_noise and ayes*snr<=nays and nays>=min_support) {
        ge.queue_node_deletion(std::abs(nv.node_id()));
        return true;
    } else {
        return false;
    }
}

void GraphContigger::solve_all_canonical(GraphEditor& ge, PairedReadsDatastore &peds, int size, bool apply){
    int total=0;
    int solved=0;
    for (auto& nv: ws.sdg.get_all_nodeviews()){
        if (nv.size()<=size and nv.is_canonical_repeat()) {
            total++;
            if (solve_canonical_repeat(ge, nv, peds)){
                solved++;
            }
        }
    }
    std::cout << solved << "/" << total << " canonical repeats solved!" << std::endl;
    if (apply) {
        ge.apply_all();
        ws.sdg.join_all_unitigs();
    }
}

void GraphContigger::clip_all_tips(GraphEditor &ge, PairedReadsDatastore &peds, int size, bool apply) {
    int total=0;
    int solved=0;
    for (auto &nv: ws.sdg.get_all_nodeviews()){
        if (nv.size()<=size and nv.is_tip()){
            total++;
            if (clip_tip(ge, nv, peds)){
                solved++;
            }
        }
    }
    std::cout << solved << "/" << total << " tips solved!" << std::endl;
    if (apply) {
        ge.apply_all();
        ws.sdg.join_all_unitigs();
    }
}

void GraphContigger::pop_all_error_bubbles(GraphEditor &ge, PairedReadsDatastore &peds, int size, bool apply){
    int total=0;
    int solved=0;
    for (auto &nv: ws.sdg.get_all_nodeviews()){
        if (nv.size()<=size and nv.is_bubble_side() and std::abs(nv.node_id())<std::abs(nv.parallels()[0].node_id())) {
            total++;
            if (pop_error_bubbble(ge, nv, nv.parallels()[0], peds)){
                solved++;
            }
        }
    }
    std::cout << solved << "/" << total << " canonical repeats solved!" << std::endl;
    if (apply) {
        ge.apply_all();
        ws.sdg.join_all_unitigs();
    }
}

bool GraphContigger::pop_error_bubbble(GraphEditor& ge, NodeView &nv1, NodeView &nv2, PairedReadsDatastore& peds, int min_support, int max_noise, int snr, bool verbose) {
//    auto nv1(_nv1);
//    auto nv2(_nv2);
    if (!nv1.is_bubble_side()) return false;

    if (nv1.parallels()[0].node_id() != nv2.node_id() or std::abs(nv1.node_id())>std::abs(nv2.node_id())) return false;
    int v1=0;
    int v2=0;
    int vo=0;
    for (const auto&p: peds.mapper.all_paths_fw(nv1.prev()[0].node().node_id(), false)) {
        if (p[0] == nv1.node_id()){
            v1++;
        } else if (p[0] == nv2.node_id()) {
            v2++;
        } else{
            vo++;
        }
    }
    for (const auto &p: peds.mapper.all_paths_fw(-nv1.next()[0].node().node_id(), false)) {
        if (p[0] == -nv1.node_id()) {
            v1++;
        } else if (p[0] == -nv2.node_id()) {
            v2++;
        } else {
            vo++;
        }
    }
    if (v1>min_support and v2+vo<max_noise and v1>=snr*(v2+vo)) {
        ge.queue_node_deletion(nv2.node_id());
        return true;
    }
    if (v2>min_support and v1+vo<max_noise and v2>=snr*(v1+vo)) {
        ge.queue_node_deletion(nv1.node_id());
        return true;
    }
    return false;
}

void GraphContigger::solve_all_tangles(WorkSpace &ws, PairedReadsDatastore &peds, int fsize, int fminkci, int fmaxkci, bool apply){
    uint64_t total=0;
    auto s=Strider(ws);
    s.add_datastore(peds);
    auto ge=GraphEditor(ws);

    std::unordered_map<std::string, uint64_t> solved;
    std::unordered_map<std::string, uint64_t> unsolved;

    for (auto &t: ws.sdg.get_all_tangleviews(fsize, fminkci, fmaxkci)) {
        auto tc = t.classify_tangle();
        if (tc == "tip"){
            if (solve_tip(t, s, ge)){
                solved[tc]++;
            } else {
                unsolved[tc]++;
            }
        } else if (tc=="bubble") {
            if (solve_bubble(t, s, ge)) {
                solved[tc]++;
            } else {
                unsolved[tc]++;
            }
        } else if (tc.find("repeat") != std::string::npos) {
            if (solve_repeat(t, s, ge)) {
                solved[tc]++;
            } else {
                unsolved[tc]++;
            }
        } else {
            if (solve_unclassified(t, s, ge)){
                solved[tc]++;
            } else {
                unsolved[tc]++;
            }
        }
        total++;
    }
    std::cout << "Total tangles: " << total << std::endl;
    std::cout << "Solved: ";
    for (const auto &k: solved) {
        std::cout << k.first << "," << k.second;
    }
    std::cout << std::endl;
    std::cout << "Unsolved: ";
    for (const auto &k: unsolved) {
        std::cout << k.first << "," << k.second;
    }
    std::cout << std::endl;

    if (apply){
        ge.apply_all();
        ws.sdg.join_all_unitigs();
    }
}

std::vector<std::string> GraphContigger::contig_reduction_to_unique_kmers(int min_cov, int max_cov, uint32_t max_run_size){
    std::vector<std::string> seqs;
    for (const auto& nv: ws.sdg.get_all_nodeviews()){
        auto c = nv.kmer_coverage("pek31", "pe");
        int i=0;
        while(i<c.size()){

            while(i<c.size() and (c[i]<min_cov or c[i]>max_cov)){
                i++;
            }
            if (i==c.size()) break;
            auto si=i;
            uint32_t run_size=0;
            while(i<c.size() and c[i]>=min_cov and c[i]<=max_cov and run_size<=max_run_size){
                run_size++;
                i++;
            }
            seqs.push_back(nv.sequence().substr(si, i-si+30));
        }
    }
    return seqs;
}

std::map<uint64_t, std::vector<sgNodeID_t >> GraphContigger::group_nodes(PairedReadsDatastore peds){
    std::cout << "Starting group nodes" << std::endl;
    std::map<uint64_t, std::vector<sgNodeID_t >> groups;
    std::unordered_set<sgNodeID_t> non_visited;
    std::unordered_set<sgNodeID_t> visited;
    std::cout << "Creating set" << std::endl;
    for (const auto &nv: ws.sdg.get_all_nodeviews()){
        non_visited.emplace(abs(nv.node_id()));
    }

    std::cout << "Starting while" << std::endl;

    while (non_visited.size()>0){

        // equivalent to set.pop()
        auto node = *non_visited.begin();
        non_visited.erase(non_visited.begin());

        auto central_node = node;
        visited.emplace(abs(node));
        auto next_node = peds.mapper.get_node_inmediate_neighbours(node);
        std::vector<sgNodeID_t > group;
        group.push_back(node);

        while (next_node != 0 and std::find(visited.begin(), visited.end(), abs(next_node)) == visited.end()){

            visited.emplace(abs(next_node));
            non_visited.erase(abs(next_node));
            group.push_back(node);
            node = next_node;
            next_node = peds.mapper.get_node_inmediate_neighbours(node);
        }
        std::cout << "Grupo: " << central_node << ": ";
        for (const auto &g: group){
            std::cout << g << ", ";
        }
        std::cout<<std::endl;

        groups[central_node] = group;
    }
    return groups;
}
