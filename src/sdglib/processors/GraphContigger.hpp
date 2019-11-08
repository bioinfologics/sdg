//
// Created by Bernardo Clavijo (EI) on 08/11/2019.
//

#pragma once
#include <sdglib/workspace/WorkSpace.hpp>


class GraphContigger {
public:
    GraphContigger (WorkSpace &_ws):ws(_ws){};
    void solve_canonical_repeats_with_single_paths(int min_support=6, int max_noise=5, float snr=10);
private:
    WorkSpace &ws;
};
