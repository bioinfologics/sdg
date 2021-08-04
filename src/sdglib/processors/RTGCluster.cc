//
// Created by Bernardo Clavijo (EI) on 04/08/2021.
//

#include "RTGCluster.hpp"

RTGCluster::RTGCluster(const ReadThreadsGraph &_rtg, int _p, int _q, int _min_node_happiness):
    rtg(_rtg), p(_p), q(-q), min_node_happiness(_min_node_happiness) {

}