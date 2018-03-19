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
    std::unordered_map<uint64_t, graphPosition> kmer_to_graphposition;
    uint k;
public:
    uniqueKmerIndex(uint k, uint memlimit = 4) : k(k) {}

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

    std::unordered_map<uint64_t, graphPosition>::const_iterator find(uint64_t hash) {
        return kmer_to_graphposition.find(hash);
    };

    std::unordered_map<uint64_t, graphPosition>::const_iterator end() {
        return kmer_to_graphposition.cend();
    };
};


#endif //BSG_UNIQUEKMERINDEX_HPP
