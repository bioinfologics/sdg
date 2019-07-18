//
// Created by Bernardo Clavijo (EI) on 2019-07-17.
//

#pragma once

#include <sdglib/graph/DistanceGraph.hpp>
#include <sdglib/workspace/WorkSpace.hpp>

/**
 * @brief Creates a new DistanceGraph with linkage information from the workspace over the nodes of an input DistanceGraph
 *
 */
class LinkageMaker {
public:
    explicit LinkageMaker(const DistanceGraph &_dg):dg(_dg){deselect_all();};
    explicit LinkageMaker(const SequenceDistanceGraph &_dg):LinkageMaker(static_cast<const DistanceGraph &>(_dg)){};

    void deselect_all();
    void select_all();
    void report_selection();
    void select_by_size(uint64_t min_size,uint64_t max_size=0);


    DistanceGraph make_topology_linkage(int radius);

    DistanceGraph make_paired_linkage(int min_reads);

    DistanceGraph make_paired_linkage_pe(int min_reads);

    //Linkage creation methods (work on selected nodes)
    std::map<std::pair<sgNodeID_t, sgNodeID_t>, uint64_t> shared_read_paths(int min_shared, std::vector<size_t> libraries, bool r1rev, bool r2rev);

    DistanceGraph make_tag_linkage(int min_tags, bool use_kmer_paths=false);

    //supporting methods
    std::vector<Link> mappings_to_multilinkage(const std::vector<LongReadMapping> &lorm_mappings, uint32_t read_size, int32_t unmapped_end=1000);

    DistanceGraph make_longRead_multilinkage(const LongReadsMapper &lorm, bool real_read_size=true, int32_t unmapped_end=1000);

    DistanceGraph make_paired10x_multilinkage(const PairedReadsMapper &prm, const LinkedReadsMapper &lirm, float min_tnscore=0.2, bool fr=false, uint64_t read_offset=0);


    const DistanceGraph &dg;
    std::vector<bool> selected_nodes;
};

