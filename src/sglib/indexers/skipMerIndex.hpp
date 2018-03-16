//
// Created by Luis Yanes (EI) on 16/03/2018.
//

#ifndef BSG_SKIPMERINDEX_HPP
#define BSG_SKIPMERINDEX_HPP


#include <sglib/SequenceGraph.h>
#include <sglib/PairedReadMapper.h>
#include <sglib/types/KmerTypes.hpp>

class SkipMerIndex {
    SequenceGraph & sg;
    std::unordered_map<uint64_t, graphPosition> kmer_to_graphposition;
    uint64_t memlimit;
    uint k;
    uint max_coverage;
public:
    SkipMerIndex(SequenceGraph &sg, uint k, uint max_coverage) :
            sg(sg), k(k), max_coverage(max_coverage)
    {
            const std::string output_prefix("./");

            SMR<KmerIDX,
            SkipMerIDXFactory<FastaRecord>,
            GraphNodeReader<FastaRecord>,
            FastaRecord,
            GraphNodeReaderParams,
            KMerIDXFactoryParams> kmerIDX_SMR({1, sg}, {k}, {memlimit, 0, max_coverage, output_prefix});

            // Get the unique_kmers from the graph into a map
            std::cout << "Indexing graph... " << std::endl;
            kmer_to_graphposition.clear();
            std::unordered_set<int32_t> seen_contigs;
            for (auto &kidx :kmerIDX_SMR.process_from_memory()) {
                kmer_to_graphposition[kidx.kmer]={kidx.contigID,kidx.pos};
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


#endif //BSG_SKIPMERINDEX_HPP
