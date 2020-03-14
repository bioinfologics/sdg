//
// Created by Bernardo Clavijo (EI) on 08/11/2019.
//

#pragma once
#include <sdglib/workspace/WorkSpace.hpp>


class GraphContigger {
public:
    GraphContigger (WorkSpace &_ws):ws(_ws){};
    void tip_clipping(int tip_size, int rounds=10);
    void remove_small_unconnected(int min_size);
    void solve_canonical_repeats_with_single_paths(const PairedReadsDatastore & prds,int min_support=6, int max_noise=5, float snr=10);
    void solve_canonical_repeats_with_paired_paths(const PairedReadsDatastore & prds,int min_support=6, int max_noise=5, float snr=10);
    void extend_to_repeats(int max_size=300);

private:
    WorkSpace &ws;
};
