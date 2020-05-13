//
// Created by Bernardo Clavijo (EI) on 01/05/2020.
//

#pragma once

#include <sdglib/workspace/WorkSpace.hpp>

class Strider {
public:

    Strider(WorkSpace & _ws, std::vector<std::string> paired_reads={}, std::vector<std::string> long_reads={});
    //walk out from a node
    SequenceDistanceGraphPath walk_out(sgNodeID_t n);

    SequenceDistanceGraphPath walk_out_in_order(sgNodeID_t n, bool use_pair=true, bool collapse_pair=true, bool verbose=false);
    WorkSpace & ws;
    std::vector<PairedReadsDatastore *> paired_datastores;
    std::vector<LongReadsDatastore *> long_datastores;
};

