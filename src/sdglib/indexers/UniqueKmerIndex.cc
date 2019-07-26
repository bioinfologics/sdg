//
// Created by Luis Yanes (EI) on 2019-06-25.
//

#include "UniqueKmerIndex.hpp"
#include <sdglib/graph/SequenceDistanceGraph.hpp>


UniqueKmerIndex::UniqueKmerIndex(const SequenceDistanceGraph &sg, const uint8_t _k) : k(_k){
    kmer_to_graphposition.clear();
    std::vector<pair> kidxv;
    uint64_t total_k { 0 };
    total_kmers_per_node = std::vector<uint64_t>(sg.nodes.size(), 0);
    for (sgNodeID_t node = 0; node < sg.nodes.size(); node++) {
        auto sgnode = sg.nodes[node];
        if (sgnode.sequence.size() >= k) {
            auto n = sgnode.sequence.size() + 1 - k;
            total_k += n;
            total_kmers_per_node[node] = n;
        }
    }
    kidxv.reserve(total_k);
    FastaRecord r;
    kmerPosFactory kcf({k});
    for (sgNodeID_t n = 1; n < sg.nodes.size(); ++n) {
        if (sg.nodes[n].sequence.size() >= k) {
            r.id = n;
            r.seq = sg.nodes[n].sequence;
            kcf.setFileRecord(r);
            kcf.next_element(kidxv);
        }
    }

    sdglib::sort(kidxv.begin(),kidxv.end(),[](const pair & a, const pair & b){return a.first<b.first;});

    auto wi=kidxv.begin();
    auto ri=kidxv.begin();
    auto nri=kidxv.begin();
    while (ri<kidxv.end()){
        while (nri!=kidxv.end() and nri->first==ri->first) ++nri;
        if (nri-ri==1) {
            *wi=*ri;
            ++wi;
        }
        ri=nri;
    }
    kidxv.resize(wi - kidxv.begin());
    std::unordered_set<sgNodeID_t > seen_contigs;
    seen_contigs.reserve(sg.nodes.size());
    unique_kmers_per_node = std::vector<uint64_t>(sg.nodes.size(), 0);
    kmer_to_graphposition.reserve(kidxv.size());
    for (auto &kidx :kidxv) {
        kmer_to_graphposition[kidx.first] = { kidx.second.node, kidx.second.pos };
        unique_kmers_per_node[std::abs(kidx.second.node)] += 1;
        seen_contigs.insert(std::abs(kidx.second.node));
    }
}

Unique63merIndex::Unique63merIndex(const SequenceDistanceGraph &sg, const uint8_t _k) : k(_k){
    kmer_to_graphposition.clear();
    std::vector<pair> kidxv;
    uint64_t total_k { 0 };
    total_kmers_per_node = std::vector<uint64_t>(sg.nodes.size(), 0);
    for (sgNodeID_t node = 0; node < sg.nodes.size(); node++) {
        auto sgnode = sg.nodes[node];
        if (sgnode.sequence.size() >= k) {
            auto n = sgnode.sequence.size() + 1 - k;
            total_k += n;
            total_kmers_per_node[node] = n;
        }
    }
    kidxv.reserve(total_k);
    FastaRecord r;
    kmerPosFactory128 kcf({k});
    for (sgNodeID_t n = 1; n < sg.nodes.size(); ++n) {
        if (sg.nodes[n].sequence.size() >= k) {
            r.id = n;
            r.seq = sg.nodes[n].sequence;
            kcf.setFileRecord(r);
            kcf.next_element(kidxv);
        }
    }

    sdglib::sort(kidxv.begin(),kidxv.end(),[](const pair & a, const pair & b){return a.first<b.first;});

    auto wi=kidxv.begin();
    auto ri=kidxv.begin();
    auto nri=kidxv.begin();
    while (ri<kidxv.end()){
        while (nri!=kidxv.end() and nri->first==ri->first) ++nri;
        if (nri-ri==1) {
            *wi=*ri;
            ++wi;
        }
        ri=nri;
    }
    kidxv.resize(wi - kidxv.begin());

    std::unordered_set<sgNodeID_t > seen_contigs;
    seen_contigs.reserve(sg.nodes.size());
    unique_kmers_per_node = std::vector<uint64_t>(sg.nodes.size(), 0);
    kmer_to_graphposition.reserve(kidxv.size());
    for (auto &kidx :kidxv) {
        kmer_to_graphposition[kidx.first] = { kidx.second.node, kidx.second.pos };
        unique_kmers_per_node[std::abs(kidx.second.node)] += 1;
        seen_contigs.insert(std::abs(kidx.second.node));
    }
}