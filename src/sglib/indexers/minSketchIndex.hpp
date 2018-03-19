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
    using Map = std::unordered_map<uint64_t, std::vector<graphStrandPos>>;
    using const_iterator = Map::const_iterator;
    using nc_pair = std::pair<uint64_t , std::vector<graphStrandPos>>;

    Map kmer_to_graphposition;
    uint8_t w = 0;
    uint8_t k = 0;
public:
    minSketchIndex(const SequenceGraph &sg, uint8_t k, uint8_t w) : k(k), w(w) {
        generate_index(sg);
    }
    minSketchIndex(uint8_t k, uint8_t w) : k(k), w(w) {
    }

    void generate_index(const SequenceGraph &sg) {
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

    const Map& getMap() {
        return kmer_to_graphposition;
    };

    const_iterator find(uint64_t hash) {
        return kmer_to_graphposition.find(hash);
    };

    const_iterator end() {
        return kmer_to_graphposition.cend();
    };

    void write_to_file(std::string &filename) {
        sglib::OutputLog() << "Dumping index" << std::endl;
        std::ofstream outf(filename, std::ios_base::binary);
        auto mapSize(kmer_to_graphposition.size());
        outf.write(reinterpret_cast<const char *>(&k), sizeof(k));
        outf.write(reinterpret_cast<const char *>(&w), sizeof(w));
        outf.write(reinterpret_cast<const char *>(&mapSize), sizeof(mapSize));
        for (const auto &skm: kmer_to_graphposition){
            outf.write(reinterpret_cast<const char *>(&skm.first), sizeof(skm.first));
            auto vSize(skm.second.size());
            outf.write(reinterpret_cast<const char *>(&vSize), sizeof(vSize));
            for(const auto &v: skm.second) {
                outf.write(reinterpret_cast<const char *>(&v), sizeof(v));
            }
        }
        sglib::OutputLog() << "Done!" << std::endl;
    }

    void load_from_file(std::string &filename) {
        sglib::OutputLog() << "Loading index" << std::endl;
        std::ifstream outf(filename, std::ios_base::binary);
        nc_pair skm;
        auto mapSize(kmer_to_graphposition.size());
        outf.read(reinterpret_cast<char *>(&k), sizeof(k));
        outf.read(reinterpret_cast<char *>(&w), sizeof(w));
        outf.read(reinterpret_cast<char *>(&mapSize), sizeof(mapSize));
        for (decltype(kmer_to_graphposition.size()) n = 0; n < mapSize; ++n){
            outf.read(reinterpret_cast<char *>(&skm.first), sizeof(skm.first));
            decltype(skm.second.size()) vSize;
            outf.read(reinterpret_cast<char *>(&vSize), sizeof(vSize));
            skm.second.resize(vSize);
            for(decltype(vSize) v = 0; v < vSize; ++v) {
                outf.read(reinterpret_cast<char *>(&skm.second[v]), sizeof(skm.second[0]));
            }
        }
        sglib::OutputLog() << "Done!" << std::endl;
    }

};


#endif //BSG_MINSKETCHINDEX_HPP
