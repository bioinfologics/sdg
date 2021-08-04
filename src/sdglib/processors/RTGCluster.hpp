//
// Created by Bernardo Clavijo (EI) on 04/08/2021.
//

#pragma once
#include <sdglib/graph/ReadThreadsGraph.hpp>

class RTGCluster {
public:
    RTGCluster(const ReadThreadsGraph & _rtg, int _p, int _q, int _min_node_happiness);

private:
    const ReadThreadsGraph & rtg;
    int p,q;
    float min_node_happiness;

};


