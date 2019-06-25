//
// Created by Luis Yanes (EI) on 2019-06-25.
//

#include "UniqueKmerIndex.hpp"
#include <sdglib/graph/SequenceDistanceGraph.hpp>

void Unique63merIndex::generate_index(const SequenceDistanceGraph &sg, bool verbose) {
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
    if (verbose) sdglib::OutputLog(sdglib::INFO)<<kidxv.size()<<" kmers in total"<<std::endl;
    if (verbose) sdglib::OutputLog(sdglib::INFO) << "  Sorting..."<<std::endl;

    sdglib::sort(kidxv.begin(),kidxv.end(),[](const pair & a, const pair & b){return a.first<b.first;});

    if (verbose) sdglib::OutputLog(sdglib::INFO) << "  Merging..."<<std::endl;
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

    if (verbose) sdglib::OutputLog(sdglib::INFO) << kidxv.size() << " unique kmers in index, creating map" << std::endl;
    std::unordered_set<sgNodeID_t > seen_contigs;
    seen_contigs.reserve(sg.nodes.size());
    unique_kmers_per_node = std::vector<uint64_t>(sg.nodes.size(), 0);
    kmer_to_graphposition.reserve(kidxv.size());
    for (auto &kidx :kidxv) {
        kmer_to_graphposition[kidx.first] = { kidx.second.node, kidx.second.pos };
        unique_kmers_per_node[std::abs(kidx.second.node)] += 1;
        seen_contigs.insert(std::abs(kidx.second.node));
    }
    if (verbose) sdglib::OutputLog(sdglib::INFO) << seen_contigs.size() << " nodes with indexed kmers" <<std::endl;
}

void UniqueKmerIndex::generate_index(const SequenceDistanceGraph &sg, bool verbose) {
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
    if (verbose) sdglib::OutputLog(sdglib::INFO)<<kidxv.size()<<" kmers in total"<<std::endl;
    if (verbose) sdglib::OutputLog(sdglib::INFO) << "  Sorting..."<<std::endl;
#ifdef _OPENMP
    __gnu_parallel::sort(kidxv.begin(),kidxv.end(),[](const pair & a, const pair & b){return a.first<b.first;});
#else
    std::sort(kidxv.begin(),kidxv.end(),[](const pair & a, const pair & b){return a.first<b.first;});
#endif

    if (verbose) sdglib::OutputLog(sdglib::INFO) << "  Merging..."<<std::endl;
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
    if (verbose) sdglib::OutputLog(sdglib::INFO) << kidxv.size() << " unique kmers in index, creating map" << std::endl;
    std::unordered_set<sgNodeID_t > seen_contigs;
    seen_contigs.reserve(sg.nodes.size());
    unique_kmers_per_node = std::vector<uint64_t>(sg.nodes.size(), 0);
    kmer_to_graphposition.reserve(kidxv.size());
    for (auto &kidx :kidxv) {
        kmer_to_graphposition[kidx.first] = { kidx.second.node, kidx.second.pos };
        unique_kmers_per_node[std::abs(kidx.second.node)] += 1;
        seen_contigs.insert(std::abs(kidx.second.node));
    }
    if (verbose) sdglib::OutputLog(sdglib::INFO) << seen_contigs.size() << " nodes with indexed kmers" <<std::endl;
}
