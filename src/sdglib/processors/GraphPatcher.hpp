//
// Created by Bernardo Clavijo (EI) on 20/05/2020.
//

#pragma once

#include <sdglib/workspace/WorkSpace.hpp>

class GraphPatcher {
public:
    GraphPatcher(WorkSpace & _ws):ws(_ws){};
    uint64_t find_tips_to_reconnect(int min_paths); //populates reconnecition_groups with 1-1 tips that connect to each other only
    void create_patch(std::vector<sgNodeID_t> reconnection_group);//takes all reads from the reconnection group, creates a local assembly, gets the lines from it.
    uint64_t collapse_reconnection_groups();
    uint64_t patch_reconnection_groups();
    std::vector<std::vector<sgNodeID_t>> reconnection_groups;
    WorkSpace & ws;
};

