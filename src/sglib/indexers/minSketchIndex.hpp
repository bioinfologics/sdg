//
// Created by Luis Yanes (EI) on 15/03/2018.
//

#ifndef BSG_MINSKETCHINDEX_HPP
#define BSG_MINSKETCHINDEX_HPP

#include <sglib/factories/StrandedMinSketchFactory.h>
#include <sglib/readers/FileReader.h>
#include <sglib/readers/SequenceGraphReader.h>
#include <cmath>
#include <sglib/types/KmerTypes.hpp>

class minSketchIndex {
    std::unordered_map<uint64_t, std::vector<graphStrandPos>> kmer_to_graphposition;
public:
    minSketchIndex(SequenceGraph &sg, uint k, uint w) {
        StrandedMinimiserSketchFactory kf(k, w);
        std::set<MinPosIDX> sketch;
        GraphNodeReader<FastaRecord> gnr({0,sg});
        FastaRecord node;
        while (gnr.next_record(node)) {
            kf.getMinSketch(node.seq, sketch);
            for (const auto &sk : sketch) {
                kmer_to_graphposition[sk.hash].emplace_back(node.id * (std::signbit(sk.pos)?-1:1), std::abs(sk.pos));
            }
            sketch.clear();
        }
    }

    std::unordered_map<uint64_t, std::vector<graphStrandPos>>::const_iterator find(uint64_t hash) {
        return kmer_to_graphposition.find(hash);
    };

    std::unordered_map<uint64_t, std::vector<graphStrandPos>>::const_iterator end() {
        return kmer_to_graphposition.cend();
    };
};


#endif //BSG_MINSKETCHINDEX_HPP
