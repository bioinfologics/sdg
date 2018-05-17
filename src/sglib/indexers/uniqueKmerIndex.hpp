//
// Created by Luis Yanes (EI) on 15/03/2018.
//

#ifndef BSG_UNIQUEKMERINDEX_HPP
#define BSG_UNIQUEKMERINDEX_HPP


#include <cstdint>
#include <unordered_set>
#include <iostream>
#include <sglib/factories/KMerIDXFactory.h>
#include <sglib/readers/SequenceGraphReader.h>
#include <sglib/readers/FileReader.h>
#include <sglib/SMR.h>
#include <sglib/PairedReadMapper.h>
#include <sglib/datastores/LinkedReadsDatastore.hpp>
#include <sglib/types/KmerTypes.hpp>

class uniqueKmerIndex {
    using Map = std::unordered_map<uint64_t, graphStrandPos>;
    using const_iterator = std::unordered_map<uint64_t, graphStrandPos>::const_iterator;
    using pair = std::pair<uint64_t, graphStrandPos>;

    Map kmer_to_graphposition;
    uint k;
    std::vector<uint64_t> unique_kmers_per_node;
    std::vector<uint64_t> total_kmers_per_node;

public:
    explicit uniqueKmerIndex(uint k) : k(k) {}

    uniqueKmerIndex(const SequenceGraph &sg, uint k, uint64_t memlimit = 4) :
            k(k)
    {
        generate_index(sg, k, memlimit);
    }

    void generate_index(const SequenceGraph &sg, uint8_t k, uint64_t memlimit = 4) {
        const std::string output_prefix("./");

        this->k = k;

        SMR<KmerIDX,
        kmerIDXFactory<FastaRecord>,
        GraphNodeReader<FastaRecord>,
        FastaRecord,
        GraphNodeReaderParams,
        KMerIDXFactoryParams> kmerIDX_SMR({1, sg}, {k}, {memlimit*GB, 0, 1, output_prefix});

        // Get the unique_kmers from the graph into a map
        std::cout << "Indexing graph... " << std::endl;
        kmer_to_graphposition.clear();
        std::unordered_set<int32_t> seen_contigs;
        unique_kmers_per_node = std::vector<uint64_t>(sg.nodes.size(), 0);
        total_kmers_per_node = std::vector<uint64_t>(sg.nodes.size(), 0);

        for (auto &kidx : kmerIDX_SMR.process_from_memory()) {
            kmer_to_graphposition[kidx.kmer] = {kidx.contigID, kidx.pos};
            unique_kmers_per_node[std::abs(kidx.contigID)] += 1;
            seen_contigs.insert(std::abs(kidx.contigID));
        }

        for (sgNodeID_t node = 0; node < sg.nodes.size(); node++) {
            total_kmers_per_node[node] = sg.nodes[node].sequence.size() - (k - 1);
        }
    }

    const_iterator find(const uint64_t hash) const {
        return kmer_to_graphposition.find(hash);
    };

    const_iterator end() const {
        return kmer_to_graphposition.cend();
    };

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

    void write_vector(std::ofstream &outf, std::vector<uint64_t> &vec) {
        auto vecSize(vec.size());
        outf.write(reinterpret_cast<const char *>(&vecSize), sizeof(vecSize));
        outf.write(reinterpret_cast<const char *>(vec.data()),
                   vecSize * sizeof(uint64_t));
    }

    void read_vector(std::ifstream &input, std::vector<uint64_t> &v) {
        auto vSize(v.size());
        input.read(reinterpret_cast<char *>(&vSize), sizeof(vSize));
        v.resize(vSize);
        input.read(reinterpret_cast<char *>(v.data()), vSize * sizeof(uint64_t));
        //v.insert(v.begin(), std::istream_iterator<uint64_t>(input), std::istream_iterator<uint64_t>());
    }

    void write_to_file(std::string filename) {
        sglib::OutputLog() << "Dumping index" << std::endl;
        std::ofstream outf(filename, std::ios_base::binary);
        outf.write(reinterpret_cast<const char *>(&k), sizeof(k));
        auto mapSize(kmer_to_graphposition.size());
        outf.write(reinterpret_cast<const char *>(&mapSize), sizeof(mapSize));
        for (const auto &skm: kmer_to_graphposition){
            outf.write(reinterpret_cast<const char *>(&skm.first), sizeof(skm.first));
            outf.write(reinterpret_cast<const char *>(&skm.second), sizeof(skm.second));
        }

        write_vector(outf, unique_kmers_per_node);
        write_vector(outf, total_kmers_per_node);

        sglib::OutputLog() << "Done!" << std::endl;
    }

    void read_from_file(std::string filename) {
        sglib::OutputLog() << "Loading index" << std::endl;
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

        read_vector(inf, unique_kmers_per_node);
        read_vector(inf, total_kmers_per_node);

        sglib::OutputLog() << "Done!" << std::endl;
    }

    bool operator==(const uniqueKmerIndex &other) const {
        return k == other.k && kmer_to_graphposition == other.kmer_to_graphposition && unique_kmers_per_node == other.unique_kmers_per_node && total_kmers_per_node == other.total_kmers_per_node;
    }
};


#endif //BSG_UNIQUEKMERINDEX_HPP
