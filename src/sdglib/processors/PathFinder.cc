//
// Created by Bernardo Clavijo (EI) on 07/04/2020.
//

#include "PathFinder.hpp"

void PathFinder::index_paths() {
    StringKMerFactoryNC skf(k);
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