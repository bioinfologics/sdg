//
// Created by Ben Ward (EI) on 07/02/2018.
//

#include <tuple>
#include <stack>
#include <sglib/processors/UniqueKmerIndex.h>
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
        kmer_to_node_map[kidx.kmer] = graphPosition{kidx.contigID, kidx.pos};
        unique_kmers_per_node[std::abs(kidx.contigID)] += 1;
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
    graphPosition p;
    if (exists) {
        p = nk->second;
    }
    return std::make_tuple(exists, p);
}

bool UniqueKmerIndex::is_unmappable(sgNodeID_t id) const {
    return 0 == unique_kmers_per_node[std::abs(id)];
}

bool UniqueKmerIndex::traverse_dark_nodes(const sgNodeID_t start, const sgNodeID_t goal) {
    // Create a stack with the nodes and the path length
    struct visitor {
        sgNodeID_t node;
        unsigned long dist;
        unsigned long path_length;
        visitor(sgNodeID_t n, uint d, uint p) : node(n), dist(d), path_length(p) {}
        bool operator<(const visitor& rhs) const {
            return std::tie(dist, path_length) < std::tie(rhs.dist, rhs.path_length);
        }
        bool operator<(const visitor& rhs) {
            return std::tie(dist, path_length) < std::tie(rhs.dist, rhs.path_length);
        }
        bool operator>(const visitor& rhs) const {
            return std::tie(dist, path_length) > std::tie(rhs.dist, rhs.path_length);
        }
        bool under_edge_limit(uint limit) const {
            return path_length < limit or limit == 0;
        }
        bool under_size_limit(uint limit) const {
            return dist < limit or limit == 0;
        }
        bool unexplored(const std::set<sgNodeID_t>& exp) const {
            return exp.find(node) == exp.end();
        }
    };


    std::stack<visitor> to_visit;
    to_visit.emplace(start, 0, 0);
    std::set<sgNodeID_t> explored;
    while (!to_visit.empty()) {
        // Consider the current node, assign it to activeNode.
        const auto activeNode(to_visit.top());
        to_visit.pop();
        // Check if the node has already been explored.
        // If it hasn't, add it to the set of explored nodes and continue.
        if (activeNode.unexplored(explored)) {
            explored.emplace(activeNode.node);
            // Now, find the children that are dark nodes, if any exist, place them on the to_visit column.
            for (const auto &l : sg.get_fw_links(activeNode.node)) {
                if (l.dest == goal) {

                } else if (is_unmappable(l.dest)) {
                    to_visit.emplace(l.dest, activeNode.dist + sg.nodes[std::abs(l.dest)].sequence.length(), activeNode.path_length + 1);
                }
            }
        }
    }
    //return std::vector<sgNodeID_t>(visited.begin(), visited.end());
    return true;
}