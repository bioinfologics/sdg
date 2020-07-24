//
// Created by Bernardo Clavijo (EI) on 08/11/2019.
//

#pragma once
#include <sdglib/workspace/WorkSpace.hpp>
#include <sdglib/mappers/LongReadsRecruiter.hpp>


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

private:
    WorkSpace &ws;
};
