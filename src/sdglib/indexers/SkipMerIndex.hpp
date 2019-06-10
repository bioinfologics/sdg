//
// Created by Luis Yanes (EI) on 16/03/2018.
//

#ifndef BSG_SKIPMERINDEX_HPP
#define BSG_SKIPMERINDEX_HPP

#include <sdglib/utilities/omp_safe.hpp>
#include <sdglib/graph/SequenceDistanceGraph.hpp>
#include <sdglib/types/KmerTypes.hpp>
#include <sdglib/readers/FileReader.hpp>
#include <sdglib/readers/SequenceGraphReader.hpp>
#include <cmath>
#include <sdglib/factories/SkipMerFactory.hpp>
#include <sdglib/utilities/omp_safe.hpp>

class SkipMerIndex {

    using Map = std::unordered_map<uint64_t, graphStrandPos>;
    using const_iterator = std::unordered_map<uint64_t, graphStrandPos>::const_iterator;
    using pair = std::pair<uint64_t, graphStrandPos>;

    Map kmer_to_graphposition;
    uint8_t k=0;
    uint8_t m=0;
    uint8_t n=0;
    std::vector<uint64_t> unique_kmers_per_node;
    std::vector<uint64_t> total_kmers_per_node;

public:
    SkipMerIndex(){}
    SkipMerIndex(uint8_t k, uint8_t m, uint8_t n) :
            k(k), m(m), n(n)
    {
    }

    SkipMerIndex(const SequenceGraph& sg, uint8_t k, uint8_t m, uint8_t n) :
            k(k), m(m), n(n)
    {
        generate_index(sg);
    }

    void generate_index(const SequenceGraph & sg, bool verbose = true) {
        sdglib::OutputLog() << "Generating index for " << sg.filename << std::endl;
        const std::string output_prefix("./");

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
        std::vector<std::pair<uint64_t,graphStrandPos>> kidxv;
        FastaRecord r;
        kidxv.reserve(total_k);
        SkipmerIndexFactory kf( k, m, n );
        for (sgNodeID_t n = 1; n < sg.nodes.size(); ++n) {
            if (sg.nodes[n].sequence.size() >= k) {
                r.id = n;
                r.seq = sg.nodes[n].sequence;
                kf.setFileRecord(r);
                kf.next_element(kidxv);
            }
        }

        if (verbose) sdglib::OutputLog(sdglib::INFO)<<kidxv.size()<<" kmers in total"<<std::endl;
        if (verbose) sdglib::OutputLog(sdglib::INFO) << "  Sorting..."<<std::endl;
        sdglib::sort(kidxv.begin(),kidxv.end(),[](const std::pair<uint64_t,graphStrandPos> & a, const std::pair<uint64_t,graphStrandPos> & b){return a.first<b.first;});

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
        sdglib::OutputLog(sdglib::INFO) << kidxv.size() << " unique kmers in index, creating map" << std::endl;
        std::unordered_set<sgNodeID_t > seen_contigs;
        seen_contigs.reserve(sg.nodes.size());
        unique_kmers_per_node = std::vector<uint64_t>(sg.nodes.size(), 0);
        for (auto &kidx :kidxv) {
            kmer_to_graphposition[kidx.first] = { kidx.second.node, kidx.second.pos };
            unique_kmers_per_node[std::abs(kidx.second.node)] += 1;
            seen_contigs.insert(std::abs(kidx.second.node));
        }
        sdglib::OutputLog(sdglib::INFO) << seen_contigs.size() << " nodes with indexed kmers" <<std::endl;

        sdglib::OutputLog() << "Done!" << std::endl;
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
        sdglib::OutputLog() << "Dumping index" << std::endl;
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
        sdglib::OutputLog() << "Done!" << std::endl;
    }

    void load_from_file(std::string &filename) {
        sdglib::OutputLog() << "Loading index" << std::endl;
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
        sdglib::OutputLog() << "Done!" << std::endl;
    }


};


#endif //BSG_SKIPMERINDEX_HPP
