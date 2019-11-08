//
// Created by Bernardo Clavijo (EI) on 08/11/2019.
//

#include "GraphContigger.hpp"
#include "GraphEditor.hpp"
#include <sdglib/views/NodeView.hpp>

void GraphContigger::solve_canonical_repeats_with_single_paths(int min_support, int max_noise, float snr) {
    GraphEditor ge(ws);
    uint64_t repeats=0;
    uint64_t solved_repeats=0;
    for (auto & nv :ws.sdg.get_all_nodeviews()){
        auto nvp=nv.prev();
        auto nvn=nv.next();
        if (nvp.size()!=2 or nvn.size()!=2) continue;
        //sdglib::OutputLog()<<"repeat on node "<<nv<<std::endl;
        ++repeats;
        std::unordered_map<std::pair<sgNodeID_t,sgNodeID_t>,uint64_t> through;
        for (auto &pid:ws.paired_reads_datastores[0].mapper.paths_in_node[nv.node_id()]){
            //sdglib::OutputLog()<<"checking path "<<pid<<std::endl;
            auto p=ws.paired_reads_datastores[0].mapper.read_paths[llabs(pid)].path;
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
    }
    sdglib::OutputLog()<<"Contigger repeat_expansion: "<<solved_repeats<<" / "<<repeats<<std::endl;
    ge.apply_all();
}