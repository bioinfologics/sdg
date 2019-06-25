//
// Created by Luis Yanes (EI) on 15/03/2018.
//

#ifndef BSG_UNIQUEKMERINDEX_HPP
#define BSG_UNIQUEKMERINDEX_HPP

#include <cstdint>
#include <unordered_map>
#include <iostream>
#include <sdglib/factories/KMerIDXFactory.hpp>
#include <sdglib/utilities/omp_safe.hpp>
#include <sdglib/types/KmerTypes.hpp>
#include <sdglib/readers/FileReader.hpp>
#include <sdglib/factories/KmerPosFactory.hpp>
#include <sdglib/utilities/io_helpers.hpp>

class SequenceDistanceGraph;
class UniqueKmerIndex {
public:
    using Map = std::unordered_map<uint64_t, graphStrandPos>;
    using pair = std::pair<uint64_t, graphStrandPos>;
    using const_iterator = std::unordered_map<uint64_t, graphStrandPos>::const_iterator;
    explicit UniqueKmerIndex(uint k = 31) : k(k) {}

    UniqueKmerIndex(const SequenceDistanceGraph &sg, uint8_t k) :
            k(k) { }

    void generate_index(const SequenceDistanceGraph &sg, bool verbose=true);

    const_iterator find(const uint64_t hash) const {
        return kmer_to_graphposition.find(hash);
    };

    const_iterator end() const {
        return kmer_to_graphposition.cend();
    };

    uint8_t get_k() const {
        return k;
    }

    std::tuple<bool, graphStrandPos> find_unique_kmer_in_graph(const uint64_t kmer) const {
        const auto nk = find(kmer);
        auto exists = end() != nk;
        graphStrandPos p;
        if (exists) {
            p = nk->second;
        }
        return std::make_tuple(exists, p);
    }

    bool is_unmappable(sgNodeID_t id) const {
        return 0 == unique_kmers_per_node[std::abs(id)];
    }

    uint64_t unique_kmers_in_node(const sgNodeID_t node) const {
        return unique_kmers_per_node[std::abs(node)];
    }

    uint64_t total_kmers_in_node(const sgNodeID_t node) const {
        return total_kmers_per_node[std::abs(node)];
    }

    void write_to_file(std::string filename) {
        sdglib::OutputLog() << "Dumping index" << std::endl;
        std::ofstream outf(filename, std::ios_base::binary);
        outf.write(reinterpret_cast<const char *>(&k), sizeof(k));
        auto mapSize(kmer_to_graphposition.size());
        outf.write(reinterpret_cast<const char *>(&mapSize), sizeof(mapSize));
        for (const auto &skm: kmer_to_graphposition){
            outf.write(reinterpret_cast<const char *>(&skm.first), sizeof(skm.first));
            outf.write(reinterpret_cast<const char *>(&skm.second), sizeof(skm.second));
        }

        sdglib::write_flat_vector(outf, unique_kmers_per_node);
        sdglib::write_flat_vector(outf, total_kmers_per_node);

        sdglib::OutputLog() << "Done!" << std::endl;
    }

    void read_from_file(std::string filename) {
        sdglib::OutputLog() << "Loading index" << std::endl;
        std::ifstream inf(filename, std::ios_base::binary);
        pair skm{0, graphStrandPos(0, 0)};
        inf.read(reinterpret_cast<char *>(&k), sizeof(k));
        auto mapSize(kmer_to_graphposition.size());
        inf.read(reinterpret_cast<char *>(&mapSize), sizeof(mapSize));
        kmer_to_graphposition.reserve(mapSize);
        for (decltype(kmer_to_graphposition.size()) n = 0; n < mapSize; ++n){
            inf.read(reinterpret_cast<char *>(&skm.first), sizeof(skm.first));
            inf.read(reinterpret_cast<char *>(&skm.second), sizeof(skm.second));
            kmer_to_graphposition.emplace(skm.first, skm.second);
        }

        sdglib::read_flat_vector(inf, unique_kmers_per_node);
        sdglib::read_flat_vector(inf, total_kmers_per_node);

        sdglib::OutputLog() << "Done!" << std::endl;
    }

    std::string kmer_to_string(uint64_t kmer, int K) {
        static char nucleotides [4] = {'A', 'C', 'G', 'T'};
        std::string s;
        for (int i = 0; i < K; ++i) {
            s += nucleotides[kmer % 4];
            kmer = kmer / 4;
        }
        return s;
    }

    void write_kmers_to_file(std::ofstream& out) {
        for (const auto& it : kmer_to_graphposition) {
            out << kmer_to_string(it.first, get_k()) << std::endl;
        }
        out.flush();
    }

    bool operator==(const UniqueKmerIndex &other) const {
        return k == other.k && kmer_to_graphposition == other.kmer_to_graphposition && unique_kmers_per_node == other.unique_kmers_per_node && total_kmers_per_node == other.total_kmers_per_node;
    }

    const Map& getMap() const {return kmer_to_graphposition; }
private:
    Map kmer_to_graphposition;
    uint8_t k = 31;
    std::vector<uint64_t> unique_kmers_per_node;
    std::vector<uint64_t> total_kmers_per_node;

};

class Unique63merIndex {
public:
    using Map = std::unordered_map<__uint128_t, graphStrandPos, int128_hash>;
    using pair = std::pair<__uint128_t, graphStrandPos>;
    using const_iterator = Map::const_iterator;
    explicit Unique63merIndex() : k(63) {}

    Unique63merIndex(const SequenceDistanceGraph &sg) :
            k(63) { }

    void generate_index(const SequenceDistanceGraph &sg, bool verbose=true);

    const_iterator find(const __uint128_t hash) const {
        return kmer_to_graphposition.find(hash);
    };

    const_iterator end() const {
        return kmer_to_graphposition.cend();
    };

    uint8_t get_k() const {
        return k;
    }

    std::tuple<bool, graphStrandPos> find_unique_kmer_in_graph(const __uint128_t kmer) const {
        const auto nk = find(kmer);
        auto exists = end() != nk;
        graphStrandPos p;
        if (exists) {
            p = nk->second;
        }
        return std::make_tuple(exists, p);
    }

    bool is_unmappable(sgNodeID_t id) const {
        return 0 == unique_kmers_per_node[std::abs(id)];
    }

    uint64_t unique_kmers_in_node(const sgNodeID_t node) const {
        return unique_kmers_per_node[std::abs(node)];
    }

    uint64_t total_kmers_in_node(const sgNodeID_t node) const {
        return total_kmers_per_node[std::abs(node)];
    }

    void write_to_file(std::string filename) {
        sdglib::OutputLog() << "Dumping index" << std::endl;
        std::ofstream outf(filename, std::ios_base::binary);
        outf.write(reinterpret_cast<const char *>(&k), sizeof(k));
        auto mapSize(kmer_to_graphposition.size());
        outf.write(reinterpret_cast<const char *>(&mapSize), sizeof(mapSize));
        for (const auto &skm: kmer_to_graphposition){
            outf.write(reinterpret_cast<const char *>(&skm.first), sizeof(skm.first));
            outf.write(reinterpret_cast<const char *>(&skm.second), sizeof(skm.second));
        }

        sdglib::write_flat_vector(outf, unique_kmers_per_node);
        sdglib::write_flat_vector(outf, total_kmers_per_node);

        sdglib::OutputLog() << "Done!" << std::endl;
    }

    void read_from_file(std::string filename) {
        sdglib::OutputLog() << "Loading index" << std::endl;
        std::ifstream inf(filename, std::ios_base::binary);
        pair skm{0, graphStrandPos(0, 0)};
        inf.read(reinterpret_cast<char *>(&k), sizeof(k));
        auto mapSize(kmer_to_graphposition.size());
        inf.read(reinterpret_cast<char *>(&mapSize), sizeof(mapSize));
        kmer_to_graphposition.reserve(mapSize);
        for (decltype(kmer_to_graphposition.size()) n = 0; n < mapSize; ++n){
            inf.read(reinterpret_cast<char *>(&skm.first), sizeof(skm.first));
            inf.read(reinterpret_cast<char *>(&skm.second), sizeof(skm.second));
            kmer_to_graphposition.emplace(skm.first, skm.second);
        }

        sdglib::read_flat_vector(inf, unique_kmers_per_node);
        sdglib::read_flat_vector(inf, total_kmers_per_node);

        sdglib::OutputLog() << "Done!" << std::endl;
    }

    std::string kmer_to_string(uint64_t kmer, int K) {
        static char nucleotides [4] = {'A', 'C', 'G', 'T'};
        std::string s;
        for (int i = 0; i < K; ++i) {
            s += nucleotides[kmer % 4];
            kmer = kmer / 4;
        }
        return s;
    }

    void write_kmers_to_file(std::ofstream& out) {
        for (const auto& it : kmer_to_graphposition) {
            out << kmer_to_string(it.first, get_k()) << std::endl;
        }
        out.flush();
    }

    const Map& getMap() const {return kmer_to_graphposition; }

    bool operator==(const Unique63merIndex &other) const {
        return k == other.k && kmer_to_graphposition == other.kmer_to_graphposition && unique_kmers_per_node == other.unique_kmers_per_node && total_kmers_per_node == other.total_kmers_per_node;
    }

private:
    Map kmer_to_graphposition;
    uint8_t k;
    std::vector<uint64_t> unique_kmers_per_node;
    std::vector<uint64_t> total_kmers_per_node;

};


#endif //BSG_UNIQUEKMERINDEX_HPP
