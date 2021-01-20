//
// Created by Bernardo Clavijo (EI) on 08/11/2019.
//

#pragma once
#include <sdglib/workspace/WorkSpace.hpp>
#include <sdglib/mappers/LongReadsRecruiter.hpp>
#include <sdglib/processors/Strider.hpp>
#include <sdglib/processors/GraphEditor.hpp>

class GraphContigger {
public:
    GraphContigger (WorkSpace &_ws):ws(_ws){};
    void reconnect_tips(const PairedReadsDatastore & prds, int min_support=6);
    void clip_tips(int tip_size, int rounds=10);
    void pop_bubbles(const PairedReadsDatastore & prds, int bubble_size, int min_support=6, int max_noise=5, float snr=10);
    void remove_small_unconnected(int min_size);
    void solve_canonical_repeats_with_single_paths(const PairedReadsDatastore & prds,int min_support=6, int max_noise=5, float snr=10, bool join_unitigs = true, bool dry_run=false, bool verbose=false);
    void solve_canonical_repeats_with_paired_paths(const PairedReadsDatastore & prds,int min_support=6, int max_noise=5, float snr=10, bool join_unitigs = true, bool dry_run=false, bool verbose=false);
    void solve_canonical_repeats_with_long_reads(const LongReadsRecruiter & lrr, float max_side_kci=1.5, int min_support=6, int max_noise=5, float snr=10);
    void extend_to_repeats(int max_size=300);

    // graph simplification
    bool solve_canonical_repeat(GraphEditor& ge, NodeView &nv, PairedReadsDatastore& peds, int min_support=5, int max_noise=10, int snr=10, bool verbose=false);
    bool clip_tip(GraphEditor& ge, const NodeView &_nv, PairedReadsDatastore& peds, int min_support=5, int max_noise=10, int snr=10, bool verbose=false);
    bool pop_error_bubbble(GraphEditor& ge, NodeView &nv1, NodeView &nv2, PairedReadsDatastore& peds, int min_support=5, int max_noise=10, int snr=10, bool verbose=false);

    void solve_all_canonical(GraphEditor& ge, PairedReadsDatastore &peds, int size=1000, bool apply=false);
    void clip_all_tips(GraphEditor& ge, PairedReadsDatastore &peds, int size=300, bool apply=false);
    void pop_all_error_bubbles(GraphEditor& ge, PairedReadsDatastore &peds, int size=300, bool apply=false);

    // Tangle resolution
    bool solve_bubble(TangleView &t, Strider &s, GraphEditor &ge);
    bool solve_repeat(TangleView &t, Strider &s, GraphEditor &ge);
    bool solve_tip(TangleView &t, Strider &s, GraphEditor &ge);
    bool solve_unclassified(TangleView &t, Strider &s, GraphEditor &ge);
    std::vector<sgNodeID_t> end_to_end_solution(std::vector<sgNodeID_t> p, std::vector<sgNodeID_t> fnids);

    void solve_all_tangles(WorkSpace &ws, PairedReadsDatastore& peds, int fsize=220, int fminkci=-1, int fmaxkci=-1, bool apply=false);

    // Unique kmers pipeline
    std::vector<std::string> contig_reduction_to_unique_kmers(std::string kmer_counter, std::string kmer_count, int min_cov, int max_cov, uint32_t max_run_size=100);
    std::map<uint64_t, std::vector<sgNodeID_t >> group_nodes(PairedReadsDatastore peds);

private:
    WorkSpace &ws;
};
