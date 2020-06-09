//
// Created by Bernardo Clavijo (EI) on 07/04/2020.
//

#include "PathFinder.hpp"
#include <sdglib/views/NodeView.hpp>

PFScoredPath::PFScoredPath(const PathFinder &_pf, sgNodeID_t _from, sgNodeID_t _to):pf(_pf),path(_pf.ws.sdg),from(_from),to(_to){

}


void PFScoredPath::find_hits() {
    //this imposes a delta change of 63 completely hardcoded. it works but needs improvement.
    //clear old hits
    read_hitpos.clear();
    read_hitpos.resize(pf.seqs.size());

    auto seq=path.sequence();
    auto first_node_remove_bases=(pf.ovl<pf.ws.sdg.get_node_size(path.nodes[0]) ? pf.ws.sdg.get_node_size(path.nodes[0])-pf.ovl:0);
    auto last_node_remove_bases=(pf.ovl<pf.ws.sdg.get_node_size(path.nodes.back()) ? pf.ws.sdg.get_node_size(path.nodes.back())-pf.ovl:0);
    seq=seq.substr(first_node_remove_bases,seq.size()-first_node_remove_bases-last_node_remove_bases);
    //std::cout<<"path seq: "<<seq<<std::endl;
    if (seq.empty() or pf.seqs.empty()) return;
    std::vector<int64_t> last_hit(read_hitpos.size(),-1);
    std::vector<int64_t> last_delta(read_hitpos.size(),-1000000);
    StringKMerFactoryNC skf(pf.k);
    std::vector<uint64_t > pk;
    skf.create_kmers(seq,pk);
    //std::cout<<"path kmers: "<<pk.size()<<std::endl;
    //std::cout<<"scored kmers in PF: "<<pf.kmerpos.size()<<std::endl;
    for (auto &hp:read_hitpos) hp.resize(pk.size(),-1);
    //this should be done per each read k-mer independently, probably finding the mode offset first?
    for (int64_t ki=0;ki<pk.size();++ki){
        if(pf.kmerpos.count(pk[ki])==0) continue;
        for (const auto & kpos:pf.kmerpos.at(pk[ki])) {
            //std::cout<<" "<<kpos.first<<" "<<kpos.second<<" (last was "<<last_hit[kpos.first]<<" : "<<(kpos.second>last_hit[kpos.first])<<", current position points to "<<read_hitpos[kpos.first][ki]<<" : "<<(read_hitpos[kpos.first][ki]==-1)<<")"<<std::endl;
            auto delta=ki-kpos.second;
            if (read_hitpos[kpos.first][ki]==-1 and kpos.second>last_hit[kpos.first] and (last_delta[kpos.first]==-1000000 or llabs(last_delta[kpos.first] - delta )<63)){
                last_hit[kpos.first]=kpos.second;
                read_hitpos[kpos.first][ki]=kpos.second;
                last_delta[kpos.first]=delta;
            }
        }
    }
}

std::pair<uint64_t, uint64_t> PFScoredPath::score(uint64_t size) {
    std::pair<uint64_t, uint64_t> sc;
    for (auto &hp:read_hitpos){
        for (uint64_t p=0;p<hp.size() and p < size;++p){
            if (hp[p]!=-1) ++sc.first;
        }
    }
    sc.second=0;
    for (auto const &s:pf.seqs) sc.second+=s.seq.size()-pf.k+1;
    sc.second-=sc.first;
    return sc;
}

void PathFinder::index_paths(std::vector<SequenceDistanceGraphPath> _paths) {
    paths=_paths;
    StringKMerFactoryNC skf(k);
    kmerpos.clear();
    for (auto i=0;i<paths.size();++i){
        std::vector<uint64_t > pk;
        auto &p=paths[i];
        skf.create_kmers(p.sequence(),pk);
        for (auto ki=0;ki<pk.size();++ki) {
            auto &kmer=pk[ki];
            if (kmerpos.count(kmer)==0) kmerpos[kmer]={};
            kmerpos[kmer].emplace_back(i,ki);
        }
    }
}

std::vector<std::vector<uint64_t> > PathFinder::seq_to_pathpos(uint16_t path_id, std::string seq) {
    StringKMerFactoryNC skf(k);
    std::vector<uint64_t > sk;
    skf.create_kmers(seq,sk);
    std::vector<std::vector<uint64_t> > r(sk.size());
    for (auto ki=0;ki<sk.size();++ki) {
        auto &kmer = sk[ki];
        if (kmerpos.count(kmer) != 0){
            for (auto &pkp:kmerpos[kmer]){
                if (pkp.first==path_id) r[ki].emplace_back(pkp.second);
            }
        }
    }
    return r;
}

void PathFinder::load_lrseqs(DistanceGraph &dg, const LongReadsRecruiter & lrr, int ovl_extension) {
    ovl=ovl_extension;
    //get all reads joining the nodes in the dg (allows for more than 1 pass)
    for (auto &l:dg.get_nodeview(n1).next()) {
        if (l.node().node_id() != n2 or l.distance()<-200) continue;
        auto rid = l.support().id;
        int64_t fstart = -1;
        int64_t rstart = -1;
        for (auto &tn:lrr.read_threads[rid]) {
            if (tn.node == n1) fstart = std::max(tn.end - ovl_extension, 0);
            if (tn.node == -n2) rstart = std::max(tn.end - ovl_extension, 0);
            if (tn.node == n2 and fstart != -1 and fstart < tn.start) {
                //std::cout<<rid<<std::endl;std::flush(std::cout);
                seqs.emplace_back(lrr.datastore.get_read_sequence(rid).substr(fstart, tn.start + 300 - fstart),
                                  PFSEType::PFLongRead, 0, rid);
                fstart = -1;
            }
            if (tn.node == -n1 and rstart != -1 and rstart < tn.start) {
                seqs.emplace_back(
                        sdglib::str_rc(lrr.datastore.get_read_sequence(rid).substr(rstart, tn.start + 300 - rstart)),
                                  PFSEType::PFLongRead, 0, rid);
                rstart = -1;
            }
        }
    }
}

std::string PathFinder::lrseqs_as_fasta() {
    std::string r;
    for (auto &lrs:seqs) {
        r+=">rid"+std::to_string(lrs.rid)+"\n"+lrs.seq+"\n";
    }
    return r;
}

void PathFinder::index_seqs() {
    StringKMerFactoryNC skf(k);
    kmerpos.clear();
    for (auto i=0;i<seqs.size();++i){
        std::vector<uint64_t > sk;
        auto &lrs=seqs[i];
        skf.create_kmers(lrs.seq,sk);
        for (auto si=0;si<sk.size();++si) {
            auto &kmer=sk[si];
            if (kmerpos.count(kmer)==0) kmerpos[kmer]={};
            kmerpos[kmer].emplace_back(i,si);
        }
    }
}