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


/**
 * Stores a look-up table of minSketch -> ±Node,position
 * Provides a couple of helper functions such as:
 * read, write -> stores the current index to disk and loads it from disk
 * find, end -> useful for lookup
 */
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

    minSketchIndex(uint8_t k, uint8_t w) : k(k), w(w) {}

    void generate_index(const SequenceGraph &sg) {
        StrandedMinimiserSketchFactory kf(k, w);
        std::unordered_set<MinPosIDX> sketch;
        GraphNodeReader<FastaRecord> gnr({0,sg});
        FastaRecord node;
        while (gnr.next_record(node)) {
            kf.getMinSketch(node.seq, sketch);
            auto skSize(sketch.size());
            for (const auto &sk : sketch) {
                kmer_to_graphposition[sk.hash].emplace_back(node.id * (std::signbit(sk.pos)?-1:1), std::abs(sk.pos));
            }
            sketch.clear();
        }

        if (sglib::OutputLogLevel >= sglib::DEBUG) {
            sglib::OutputLog(false) << "Index stats: " << std::endl;
            sglib::OutputLog(false) << "K: " << int(k) << std::endl;
            sglib::OutputLog(false) << "W: " << int(w) << std::endl;
            std::vector<unsigned int> kmer_per_node(sg.nodes.size(), 0);
            float max_nodes_per_kmer(0);
            float avg_nodes_per_kmer(0);
            for (const auto &kmer_nodes:kmer_to_graphposition) {
                sglib::OutputLog(false) << "#[" << kmer_nodes.first << "] " << kmer_nodes.second.size();
                sglib::OutputLog(false) << " ";
                max_nodes_per_kmer = std::max(max_nodes_per_kmer, (float)kmer_nodes.second.size());
                avg_nodes_per_kmer += kmer_nodes.second.size();
                for (const auto &node:kmer_nodes.second) {
                    sglib::OutputLog(false) << node.node << " ";
                    kmer_per_node[std::abs(node.node)]++;
                }
                sglib::OutputLog(false) << std::endl;
            }
            avg_nodes_per_kmer/=kmer_to_graphposition.size();

            sglib::OutputLog(false) << std::endl;
            sglib::OutputLog(false) << "Number of kmers per node" << std::endl;
            for (sgNodeID_t node=1; node<sg.nodes.size(); ++node) {
                sglib::OutputLog(false) << node << " "
                        << sg.oldnames[node] << " "
                        << sg.nodes[node].sequence.size() - kmer_per_node[node] << " "
                        << kmer_per_node[node] << " "
                        << sg.nodes[node].sequence.size() << "\n";
            }

            sglib::OutputLog(false) << "Tot #kmers " << kmer_to_graphposition.size() << std::endl;
            sglib::OutputLog(false) << "Avg nodes x kmer " << avg_nodes_per_kmer << std::endl;
            sglib::OutputLog(false) << "Max nodes x kmer " << max_nodes_per_kmer << std::endl;
        }
        exit(0);
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
