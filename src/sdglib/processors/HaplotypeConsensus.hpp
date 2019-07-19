//
// Created by Bernardo Clavijo (EI) on 2019-02-27.
//

#ifndef BSG_HAPLOTYPECONSENSUS_HPP
#define BSG_HAPLOTYPECONSENSUS_HPP

#include <sdglib/workspace/WorkSpace.hpp>
#include <sdglib/graph/DistanceGraph.hpp>
#include <sdglib/utilities/io_helpers.hpp>
#include <iterator>

/**
 * This class creates a haplotype consensus from a backbone.
 *
 * TODO: add intermediate nodes (finds well suported nodes in this haplotype between anchors)
 * TODO: consensus sequence creates a consensus sequence for exporting
 * TODO: path description describes a path on the graph that follows this haplotype
 * TODO: create filling nodes: check where no current node fills a gap between two anchors and re-constructs sequence from primary evidence.
 * TODO: graph update: materialises all supported contiguous paths, effectively disentangling the graph
 *
 */
class HaplotypeConsensus {
public:
    HaplotypeConsensus(WorkSpace &_ws, const std::vector<std::vector<LongReadMapping>> filtered_read_mappings, const DistanceGraph &_mldg, const DistanceGraph &_ldg, const std::vector<sgNodeID_t> _backbone, const ReadPathParams &read_path_params);

    void orient_read_paths() {
#pragma omp parallel for
        for (uint32_t read = 0; read < read_cache.size(); ++read) {
            orient_read_path(read_cache[read].id);
        }
    }
    void orient_read_path(uint64_t rid);
    void build_line_path(int min_votes, int min_path_nodes);
    std::string consensus_sequence(int disconnected_distance, int min_distance);

    std::string generate_consensus(int min_votes=1, int min_path_nodes=2, int disconnected_distance=200, int min_distance=100) {
        std::cout << "Backbone nodes: " << std::endl;
        for (const auto &n: backbone) {
            if (n != 0){
                std::cout << "seq"<<std::abs(n)<<",";
            }
        }
        std::cout << std::endl;
        

        sdglib::OutputLog() << "Created read paths for " << read_cache.size() << " reads" << std::endl;
        std::cout << "Read ids: ";
        std::copy(read_cache.cbegin(), read_cache.cend(), std::ostream_iterator<ReadCacheItem>(std::cout, ", "));
        orient_read_paths();
        sdglib::OutputLog() << "Done orienting " << read_cache.size() << " read paths" << std::endl;

        build_line_path(min_votes, min_path_nodes);

        return consensus_sequence(disconnected_distance, min_distance);

    }

    void use_long_reads_from_file(std::string filename);

    void write_read_paths(std::string filename) {
        std::ofstream ofs(filename, std::ios_base::binary);
        sdglib::write_flat_vectorvector(ofs, oriented_read_paths);
    }

    void read_read_paths(std::string filename) {
        std::ifstream ifs(filename, std::ios_base::binary);
        sdglib::read_flat_vectorvector(ifs,oriented_read_paths);
    }

    bool reads_from_file;
    std::vector<std::vector<sgNodeID_t >> oriented_read_paths;

    std::map<uint64_t,std::string> read_seqs;
    std::vector<sgNodeID_t> backbone_filled_path;
    WorkSpace &ws;
    const DistanceGraph &mldg;
    const DistanceGraph &ldg;
    std::vector<ReadCacheItem> read_cache;
    const std::vector<sgNodeID_t> backbone;
};


#endif //BSG_HAPLOTYPECONSENSUS_HPP
