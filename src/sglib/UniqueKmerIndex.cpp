//
// Created by Ben Ward (EI) on 07/02/2018.
//

#include <tuple>
#include <sglib/UniqueKmerIndex.h>
#include <sglib/SMR.h>
#include <sglib/factories/KMerIDXFactory.h>
#include <sglib/readers/FileReader.h>
#include <sglib/readers/SequenceGraphReader.h>

UniqueKmerIndex::UniqueKmerIndex(const SequenceGraph& sg, const uint8_t k) : sg(sg) {

    SMR<KmerIDX, kmerIDXFactory<FastaRecord>,
    GraphNodeReader<FastaRecord>, FastaRecord,
    GraphNodeReaderParams, KMerIDXFactoryParams> kmerIDX_SMR({1, sg}, {k}, {10 * GB, 0, 1, "default"});

    sglib::OutputLog(sglib::LogLevels::INFO) << "Generating unique graph " << int(k) << "mers and building unique kmer index." << std::endl;

    unique_kmers_per_node = std::vector<uint64_t>(sg.nodes.size(), 0);
    total_kmers_per_node = std::vector<uint64_t>(sg.nodes.size(), 0);

    for (auto &kidx :kmerIDX_SMR.process_from_memory()) {
        kmer_to_node_map[kidx.kmer] = { kidx.contigID, kidx.pos };
        unique_kmers_per_node[std::abs(kidx.contigID)] += 1; // TODO: Ask bj if KmerCompression Index is useful for this, or is overkill?
    }

    for (sgNodeID_t node = 0; node < sg.nodes.size(); node++) {
        total_kmers_per_node[node] = sg.nodes[node].sequence.size() - (k - 1);
    }

    std::vector<uint64_t> uniqKmer_statistics(kmerIDX_SMR.summaryStatistics());
    sglib::OutputLog(sglib::LogLevels::INFO) << "Number of " << int(k) << "-kmers seen in assembly " << uniqKmer_statistics[0] << '.' << std::endl;
    sglib::OutputLog(sglib::LogLevels::INFO) << "Number of contigs from the assembly " << uniqKmer_statistics[2] << '.' << std::endl;

}

std::tuple<bool, graphPosition> UniqueKmerIndex::find_unique_kmer_in_graph(uint64_t kmer) const {
    auto nk = kmer_to_node_map.find(kmer);
    auto exists = kmer_to_node_map.end() != nk;
    graphPosition p{};
    if (exists) {
        p = nk->second;
    }
    return std::make_tuple(exists, p);
}