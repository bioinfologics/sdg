//
// Created by Luis Yanes (EI) on 16/03/2018.
//

#ifndef BSG_SKIPMERINDEX_HPP
#define BSG_SKIPMERINDEX_HPP


#include <sglib/SequenceGraph.h>
#include <sglib/PairedReadMapper.h>
#include <sglib/types/KmerTypes.hpp>
#include <sglib/factories/SkipMerIndexFactory.hpp>

class SkipMerIndex {
    SequenceGraph & sg;
    std::unordered_map<uint64_t, std::vector<graphStrandPos>> kmer_to_graphposition;
    uint k;
    uint m;
    uint n;
    uint max_coverage;
public:
    SkipMerIndex(SequenceGraph &sg, uint8_t k, uint8_t m, uint8_t n, uint max_coverage) :
            sg(sg), k(k), m(m), n(n), max_coverage(max_coverage)
    {
            const std::string output_prefix("./");
        SkipMerIDXFactory<FastaRecord> kf( k, m, n );
        std::vector<MinPosIDX> mers;
        GraphNodeReader<FastaRecord> gnr({0,sg});
        FastaRecord node;

        while (gnr.next_record(node)) {
            kf.setFileRecord(node);
            kf.next_element(mers);
            for (const auto &sk : mers) {
                kmer_to_graphposition[sk.hash].emplace_back(node.id * (std::signbit(sk.pos)?-1:1), std::abs(sk.pos));
            }
            mers.clear();
        }
    }

    std::unordered_map<uint64_t, std::vector<graphStrandPos>>::const_iterator find(uint64_t hash) {
        return kmer_to_graphposition.find(hash);
    };

    std::unordered_map<uint64_t, std::vector<graphStrandPos>>::const_iterator end() {
        return kmer_to_graphposition.cend();
    };
};


#endif //BSG_SKIPMERINDEX_HPP
