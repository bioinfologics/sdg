//
// Created by Bernardo Clavijo (EI) on 07/04/2020.
//

#include "PathFinder.hpp"
#include <sdglib/views/NodeView.hpp>

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
    //get all reads joining the nodes in the dg (allows for more than 1 pass)
    for (auto &l:dg.get_nodeview(n1).next()){
        if (l.node().node_id()!=n2) continue;
        auto rid=l.support().id;
        int64_t fstart=-1;
        int64_t rstart=-1;
        for (auto &tn:lrr.read_threads[rid]) {
            if (tn.node==n1) fstart=std::max(tn.end-ovl_extension,0);
            if (tn.node==-n2) rstart=std::max(tn.end-ovl_extension,0);
            if (tn.node==n2 and fstart!=-1) {
                PFSequenceEvidence(lrr.datastore.get_read_sequence(rid).substr(fstart,tn.start+300-fstart),PFSEType::PFLongRead,0,rid);
                fstart=-1;
            }
            if (tn.node==-n1 and rstart!=-1) {
                PFSequenceEvidence(strRC(lrr.datastore.get_read_sequence(rid).substr(rstart,tn.start+300-rstart)),PFSEType::PFLongRead,0,rid);
                fstart=-1;
            }
        }
    }
    //get all mappings from lrr, chop the reads between the
}