//
// Created by Bernardo Clavijo (EI) on 2019-08-22.
//

#pragma once

#include <sdglib/graph/DistanceGraph.hpp>
#include <sdglib/types/MappingTypes.hpp>

class GraphSelfAligner {

public:
    GraphSelfAligner(const DistanceGraph &_dg, int _k=31, int _max_kfreq=200):dg(_dg),k(_k),max_kfreq(_max_kfreq){};
    explicit GraphSelfAligner(const SequenceDistanceGraph &_dg, int _k=31, int _max_kfreq=200):
        GraphSelfAligner(static_cast<const DistanceGraph &>(_dg),_k,_max_kfreq){};

    void self_align();

    const DistanceGraph &dg;
    int k;
    int max_kfreq;
    std::vector<std::vector<SequenceMatch>> matches;
};