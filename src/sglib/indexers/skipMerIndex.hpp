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
    using Map = std::unordered_map<uint64_t, graphStrandPos>;
    using const_iterator = std::unordered_map<uint64_t, graphStrandPos>::const_iterator;
    using pair = std::pair<uint64_t, graphStrandPos>;

    Map kmer_to_graphposition;
    uint8_t k=0;
    uint8_t m=0;
    uint8_t n=0;
public:
    SkipMerIndex(){}
    SkipMerIndex(uint8_t k, uint8_t m, uint8_t n) :
            k(k), m(m), n(n)
    {
    }

    void generate_index(SequenceGraph & sg) {
        sglib::OutputLog() << "Generating index for " << sg.filename << std::endl;
        const std::string output_prefix("./");
        SkipMerIDXFactory<FastaRecord> kf( k, m, n );
        std::vector<MinPosIDX> mers;
        GraphNodeReader<FastaRecord> gnr({0,sg});
        FastaRecord node;

        while (gnr.next_record(node)) {
            kf.setFileRecord(node);
            kf.next_element(mers);
            for (const auto &sk : mers) {
                kmer_to_graphposition.emplace(sk.hash, graphStrandPos(node.id * (std::signbit(sk.pos)?-1:1), std::abs(sk.pos)));
            }
            mers.clear();
        }
        sglib::OutputLog() << "Done!" << std::endl;
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
        outf.write(reinterpret_cast<const char *>(&m), sizeof(m));
        outf.write(reinterpret_cast<const char *>(&n), sizeof(n));
        outf.write(reinterpret_cast<const char *>(&mapSize), sizeof(mapSize));
        for (const auto &skm: kmer_to_graphposition){
            outf.write(reinterpret_cast<const char *>(&skm.first), sizeof(skm.first));
            outf.write(reinterpret_cast<const char *>(&skm.second), sizeof(skm.second));
        }
        sglib::OutputLog() << "Done!" << std::endl;
    }

    void load_from_file(std::string &filename) {
        sglib::OutputLog() << "Loading index" << std::endl;
        std::ifstream outf(filename, std::ios_base::binary);
        pair skm{0,graphStrandPos(0,0)};
        auto mapSize(kmer_to_graphposition.size());
        kmer_to_graphposition.reserve(mapSize);
        outf.read(reinterpret_cast<char *>(&k), sizeof(k));
        outf.read(reinterpret_cast<char *>(&m), sizeof(m));
        outf.read(reinterpret_cast<char *>(&n), sizeof(n));
        outf.read(reinterpret_cast<char *>(&mapSize), sizeof(mapSize));
        for (decltype(kmer_to_graphposition.size()) n = 0; n < mapSize; ++n){
            outf.read(reinterpret_cast<char *>(&skm.first), sizeof(skm.first));
            outf.read(reinterpret_cast<char *>(&skm.second), sizeof(skm.second));
            kmer_to_graphposition.emplace(skm.first, skm.second);
        }
        sglib::OutputLog() << "Done!" << std::endl;
    }


};


#endif //BSG_SKIPMERINDEX_HPP
