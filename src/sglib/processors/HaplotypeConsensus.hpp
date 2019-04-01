//
// Created by Bernardo Clavijo (EI) on 2019-02-27.
//

#ifndef BSG_HAPLOTYPECONSENSUS_HPP
#define BSG_HAPLOTYPECONSENSUS_HPP

#include <sglib/workspace/WorkSpace.hpp>
#include <sglib/graph/LinkageDiGraph.hpp>
#include <sglib/utilities/io_helpers.hpp>

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
    HaplotypeConsensus(WorkSpace &_ws, const LinkageDiGraph &_mldg, const LinkageDiGraph &_ldg, const std::vector<sgNodeID_t> &_backbone, const ReadPathParams &read_path_params):
    ws(_ws)
    ,mldg(_mldg)
    ,ldg(_ldg)
    ,backbone(_backbone)
    {
        auto read_cache = ws.long_read_mappers[0].create_read_paths(backbone,read_path_params);
        if (!is_sorted(read_cache.begin(), read_cache.end())){
            sglib::sort(read_cache.begin(), read_cache.end());
        }
        auto max_rid = read_cache.back();
        oriented_read_paths.resize(max_rid.id + 1);

    }

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
        sglib::OutputLog() << "Created read paths for " << read_cache.size() << " reads" << std::endl;
        std::cout << "Read ids: ";
        std::copy(read_cache.cbegin(), read_cache.cend(), std::ostream_iterator<ReadCacheItem>(std::cout, ", "));
        orient_read_paths();
        sglib::OutputLog() << "Done orienting " << read_cache.size() << " read paths" << std::endl;

        build_line_path(min_votes, min_path_nodes);

        return consensus_sequence(disconnected_distance, min_distance);

    }

    void use_long_reads_from_file(std::string filename);

    void write_read_paths(std::string filename) {
        std::ofstream ofs(filename, std::ios_base::binary);
        sglib::write_flat_vectorvector(ofs, oriented_read_paths);
    }

    void read_read_paths(std::string filename) {
        std::ifstream ifs(filename, std::ios_base::binary);
        sglib::read_flat_vectorvector(ifs,oriented_read_paths);
    }

    bool reads_from_file;
    std::vector<std::vector<sgNodeID_t >> oriented_read_paths;

    std::map<uint64_t,std::string> read_seqs;
    std::vector<sgNodeID_t> backbone_filled_path;
    WorkSpace &ws;
    const LinkageDiGraph &mldg;
    const LinkageDiGraph &ldg;
    std::vector<ReadCacheItem> read_cache;
    const std::vector<sgNodeID_t> &backbone;
};


#endif //BSG_HAPLOTYPECONSENSUS_HPP
