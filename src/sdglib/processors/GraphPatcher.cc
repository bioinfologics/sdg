//
// Created by Bernardo Clavijo (EI) on 20/05/2020.
//

#include <sdglib/views/NodeView.hpp>
#include "GraphPatcher.hpp"
#include "GraphMaker.hpp"

std::string kts(uint64_t kmer, uint8_t k) {
    std::string seq;
    seq.reserve(k);
    for (int shift=(k-1)*2;shift>=0;shift-=2) {
        //std::cout<<"kmer: "
        switch ((kmer>>shift)%4){
            case 0:
                seq.push_back('A');
                break;
            case 1:
                seq.push_back('C');
                break;
            case 2:
                seq.push_back('G');
                break;
            case 3:
                seq.push_back('T');
                break;
        }
    }
    return seq;
}

uint64_t GraphPatcher::find_tips_to_reconnect(int min_paths) {
    for (auto fnv:ws.sdg.get_all_nodeviews())
        for (auto nv: {fnv,fnv.rc()}) {
            if (nv.next().empty()) {
                std::unordered_map<sgNodeID_t,uint64_t> votes;
                for (auto &fwp: ws.paired_reads_datastores[0].mapper.all_paths_fw(nv.node_id())) {
                    if (fwp[0]==0) ++votes[fwp[1]];
                    else ++votes[fwp[0]];
                }
                std::vector<sgNodeID_t> rg;
                rg.emplace_back(nv.node_id());
                for (auto &v:votes) {
                    if (v.second>min_paths) rg.emplace_back(-v.first);
                }
                if (rg.size()>1) {
                    reconnection_groups.emplace_back(rg);
                }
            }
        }
    return reconnection_groups.size();
}

uint64_t GraphPatcher::collapse_reconnection_groups() {
    std::vector<uint64_t> rgg(reconnection_groups.size());
    std::vector<std::vector<sgNodeID_t>> new_rgs;
    for (auto i=0;i<rgg.size();++i) rgg[i]=i;
    for (auto i=0;i<reconnection_groups.size();++i) {
        if (rgg[i]!=i) continue;
        std::set<int64_t> nodes_in_group(reconnection_groups[i].begin(),reconnection_groups[i].end());
        bool updated=true;
        while (updated==true) {
            updated=false;
            for (auto j = i+1; j < reconnection_groups.size(); ++j) {
                if (rgg[j]!=j) continue;
                for (auto &n:reconnection_groups[j]) {
                    if (nodes_in_group.count(n) > 0) {
                        nodes_in_group.insert(reconnection_groups[j].begin(),reconnection_groups[j].end());
                        rgg[j]=i;
                        updated=true;
                    }
                }
            }
        }
        new_rgs.emplace_back(nodes_in_group.begin(),nodes_in_group.end());
    }
    for (auto &rg:reconnection_groups){
        if (rg.size()==3){
            if (ws.sdg.get_nodeview(rg[0]).next().size()==1 and
            (ws.sdg.get_nodeview(rg[0]).next()[0].node().node_id()==rg[1] or
             ws.sdg.get_nodeview(rg[0]).next()[0].node().node_id()==rg[2]) ) rg={rg[1],rg[2]};
            else if (ws.sdg.get_nodeview(rg[1]).next().size()==1 and
                (ws.sdg.get_nodeview(rg[1]).next()[0].node().node_id()==rg[0] or
                 ws.sdg.get_nodeview(rg[1]).next()[0].node().node_id()==rg[2]) ) rg={rg[0],rg[2]};
            else if (ws.sdg.get_nodeview(rg[2]).next().size()==1 and
                (ws.sdg.get_nodeview(rg[2]).next()[0].node().node_id()==rg[0] or
                 ws.sdg.get_nodeview(rg[2]).next()[0].node().node_id()==rg[1]) ) rg={rg[0],rg[1]};
        }
    }
    std::swap(reconnection_groups,new_rgs);
    return reconnection_groups.size();
}



uint64_t GraphPatcher::patch_reconnection_groups() {
    uint64_t solved=0;
    std::set<uint64_t> to_delete;
    for (auto &rg:reconnection_groups){
        if (rg.size()==2 and ws.sdg.get_nodeview(rg[0]).next().empty() and ws.sdg.get_nodeview(rg[1]).next().empty()){
            auto ns1=ws.sdg.get_node_sequence(rg[0]);
            auto ns2=ws.sdg.get_node_sequence(-rg[1]); //node 2 is reversed, as to make it a destination
            //Find a perfect overlap between 15 and 126 bp
            uint64_t ovl;
            for (ovl=15;ovl<127;++ovl){
                if (ovl>ns1.size() or ovl>ns2.size()) continue;
                if (ns1.substr(ns1.size()-ovl)==ns2.substr(0,ovl)) break;
            }
            if (ovl<127) {
                auto nnid=ws.sdg.add_node(ns1+ns2.substr(ovl));
                for (auto li:ws.sdg.get_bw_links(rg[0])) ws.sdg.add_link(nnid,li.dest,li.dist,li.support);
                for (auto li:ws.sdg.get_bw_links(rg[1])) ws.sdg.add_link(-nnid,li.dest,li.dist,li.support);
                to_delete.insert(llabs(rg[0]));
                to_delete.insert(llabs(rg[1]));
                ++solved;
            }
        }
    }
    for (auto &n:to_delete) ws.sdg.remove_node(n);
    return solved;
}

void GraphPatcher::create_patch(std::vector<sgNodeID_t> reconnection_group) {

    //1) create a DBG wil ALL kmers from reads that have one of the last X kmers in the node and are pathed to it.

    std::unordered_set<uint64_t> goodish_kmers;
    int patch_K=21;
    std::set<uint64_t>rids;
    //First get all the reads from the nodes that have a path forward (this should really use the offsets, but they're not saved yet?)
    for (auto nid:reconnection_group) {
        auto last_node_kmer=sdglib::str_to_kmers(ws.sdg.get_node_sequence(nid).substr(ws.sdg.get_node_size(nid)-patch_K),patch_K).back().second;
        for (auto rid:ws.paired_reads_datastores[0].mapper.paths_in_node[llabs(nid)]){
            if (nid<0) rid=-rid;
            if (rids.count(llabs(rid))==0){
                auto rkmers=sdglib::str_to_kmers(ws.paired_reads_datastores[0].get_read_sequence(llabs(rid)),patch_K);
                if (std::find_if(rkmers.begin(),rkmers.end(),[last_node_kmer](auto x ){return last_node_kmer==x.second;})!=rkmers.end()){
                    for (auto &kmer:rkmers) goodish_kmers.insert(kmer.second);
                    rids.insert(llabs(rid));
                    if ((nid>0)==(rid>0) and rids.count(llabs(ws.paired_reads_datastores[0].get_read_pair(rid)))==0){
                        rids.insert(llabs(ws.paired_reads_datastores[0].get_read_pair(rid)));
                        for (auto &kmer:sdglib::str_to_kmers(ws.paired_reads_datastores[0].get_read_sequence(llabs(ws.paired_reads_datastores[0].get_read_pair(rid))),patch_K)) goodish_kmers.insert(kmer.second);
                    }
                }

            }
        }
    }
    WorkSpace tempws;
    GraphMaker gm(tempws.sdg);
    std::cout<<"creating graph"<<std::endl;
    std::vector<uint64_t> k64mers(goodish_kmers.begin(),goodish_kmers.end());
    std::sort(k64mers.begin(),k64mers.end());
    for (auto i=0;i<k64mers.size()-1;++i) {
        if (k64mers[i]>=k64mers[i+1]) {
            std::cout<<k64mers[i]<<">="<<k64mers[i+1]<<" @ "<<i<<std::endl;
            std::cout<<kts(k64mers[i],patch_K)<<">="<<kts(k64mers[i+1],patch_K)<<" @ "<<i<<std::endl;
        }
        if (kts(k64mers[i],patch_K)>=kts(k64mers[i+1],patch_K)) {
            std::cout<<k64mers[i]<<">="<<k64mers[i+1]<<" @ "<<i<<std::endl;
            std::cout<<kts(k64mers[i],patch_K)<<">="<<kts(k64mers[i+1],patch_K)<<" @ "<<i<<std::endl;
        }
    }
    gm.new_graph_from_kmerlist_64(k64mers,patch_K);

    std::cout<<"writing graph"<<std::endl;
    tempws.sdg.write_to_gfa1("tempws_sdg.gfa");
    std::cout<<"graph dumped"<<std::endl;

    //2) Path the reads (every kmer form the read is in a node, so ti should be trivial, just find the start and follow through
    //This is an ad-hoc index, and I don't care:
    std::unordered_map<uint64_t,sgNodeID_t> kmer_nodes;
    for (auto &nv:tempws.sdg.get_all_nodeviews()){
        for (auto kmer:sdglib::str_to_kmers(nv.sequence(),patch_K)) kmer_nodes[kmer.second]=(kmer.first ? nv.node_id():-nv.node_id());
    }
    std::vector<std::vector<sgNodeID_t >> read_paths;
    for(auto rid:rids){
        read_paths.emplace_back();
        auto &rp=read_paths.back();
        for (auto kmer:sdglib::str_to_kmers(ws.paired_reads_datastores[0].get_read_sequence(rid),patch_K)) {
            auto n=kmer_nodes[kmer.second];
            if (not kmer.first) n=-n;
            if (rp.empty() or rp.back()!=n) rp.emplace_back(n);
        }
    }

    for (auto &rp:read_paths){
        std::cout<<"[ ";
        for (auto n:rp) std::cout<<n<<", ";
        std::cout<<"]"<<std::endl;
    }

    //3) stride in this graph
    for (auto nid:reconnection_group) {
        auto last_node_kmer = sdglib::str_to_kmers(ws.sdg.get_node_sequence(nid).substr(ws.sdg.get_node_size(nid) - patch_K), patch_K).back();
        std::cout<<" Entering ro reconnect node "<<nid<<" lands on temp node "<<kmer_nodes[last_node_kmer.second]*(last_node_kmer.first? 1:-1)<<std::endl;
    }
}