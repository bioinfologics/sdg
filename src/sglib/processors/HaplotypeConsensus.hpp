//
// Created by Bernardo Clavijo (EI) on 2019-02-27.
//

#ifndef BSG_HAPLOTYPECONSENSUS_HPP
#define BSG_HAPLOTYPECONSENSUS_HPP

#include <sglib/workspace/WorkSpace.hpp>
#include <sglib/graph/LinkageDiGraph.hpp>

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
    HaplotypeConsensus(WorkSpace &_ws, const LinkageDiGraph &_mldg, const LinkageDiGraph &_ldg, std::vector<sgNodeID_t> _backbone):ws(_ws),mldg(_mldg),ldg(_ldg){
        backbone=_backbone;
        for (auto n:backbone){
            for (auto rin:ws.long_read_mappers[0].reads_in_node[llabs(n)]) long_reads_in_backbone.insert(rin);
        }
    }

    void orient_read_path(uint64_t rid);
    void build_line_path();
    std::string consensus_sequence();

    void use_long_reads_from_file(std::string filename);

    bool reads_from_file;
    std::vector<std::vector<sgNodeID_t >> oriented_read_paths;
    std::map<uint64_t,std::string> read_seqs;
    std::vector<sgNodeID_t> backbone;
    std::vector<sgNodeID_t> backbone_filled_path;
    std::set<uint64_t> long_reads_in_backbone;
    WorkSpace &ws;
    const LinkageDiGraph &mldg;
    const LinkageDiGraph &ldg;
};


#endif //BSG_HAPLOTYPECONSENSUS_HPP
