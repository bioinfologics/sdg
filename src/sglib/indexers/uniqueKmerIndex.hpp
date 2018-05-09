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
public:
    explicit uniqueKmerIndex(uint k, uint memlimit = 4) : k(k) {}

    uniqueKmerIndex(const SequenceGraph &sg, uint k, uint64_t memlimit = 4) :
            k(k)
    {
        generate_index(sg, k, memlimit);
    }

    void generate_index(const SequenceGraph &sg, uint k, uint64_t memlimit) {
        const std::string output_prefix("./");

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
        for (auto &kidx :kmerIDX_SMR.process_from_memory()) {
            kmer_to_graphposition[kidx.kmer]={kidx.contigID, kidx.pos};
            seen_contigs.insert((kidx.contigID>0?kidx.contigID:-kidx.contigID));
        }
    }

    const_iterator find(const uint64_t hash) {
        return kmer_to_graphposition.find(hash);
    };

    const_iterator end() {
        return kmer_to_graphposition.cend();
    };

    void write_to_file(std::string& filename) {
        sglib::OutputLog() << "Dumping index" << std::endl;
        std::ofstream outf(filename, std::ios_base::binary);
        auto mapSize(kmer_to_graphposition.size());
        outf.write(reinterpret_cast<const char *>(&k), sizeof(k));
        outf.write(reinterpret_cast<const char *>(&mapSize), sizeof(mapSize));
        for (const auto &skm: kmer_to_graphposition){
            outf.write(reinterpret_cast<const char *>(&skm.first), sizeof(skm.first));
            outf.write(reinterpret_cast<const char *>(&skm.second), sizeof(skm.second));
        }
        sglib::OutputLog() << "Done!" << std::endl;

    }

    void read_from_file(std::string& filename) {
        sglib::OutputLog() << "Loading index" << std::endl;
        std::ifstream outf(filename, std::ios_base::binary);
        pair skm{0,graphStrandPos(0,0)};
        auto mapSize(kmer_to_graphposition.size());
        kmer_to_graphposition.reserve(mapSize);
        outf.read(reinterpret_cast<char *>(&k), sizeof(k));
        outf.read(reinterpret_cast<char *>(&mapSize), sizeof(mapSize));
        for (decltype(kmer_to_graphposition.size()) n = 0; n < mapSize; ++n){
            outf.read(reinterpret_cast<char *>(&skm.first), sizeof(skm.first));
            outf.read(reinterpret_cast<char *>(&skm.second), sizeof(skm.second));
            kmer_to_graphposition.emplace(skm.first, skm.second);
        }
        sglib::OutputLog() << "Done!" << std::endl;
    }
};


#endif //BSG_UNIQUEKMERINDEX_HPP
