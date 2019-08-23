//
// Created by Bernardo Clavijo (EI) on 2019-08-22.
//

#pragma once

#include <sdglib/graph/DistanceGraph.hpp>
#include <sdglib/types/MappingTypes.hpp>

class GraphSelfAligner {

public:
    GraphSelfAligner(const DistanceGraph &_dg):dg(_dg){};
    explicit GraphSelfAligner(const SequenceDistanceGraph &_dg):GraphSelfAligner(static_cast<const DistanceGraph &>(_dg)){};

    void self_align();

    const DistanceGraph &dg;
    std::vector<std::vector<SequenceMatch>> matches;
};