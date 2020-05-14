//
// Created by Bernardo Clavijo (EI) on 01/05/2020.
//

#pragma once

#include <sdglib/workspace/WorkSpace.hpp>
#include <sdglib/mappers/LongReadsRecruiter.hpp>

class Strider {
public:

    Strider(WorkSpace & _ws);
    void add_datastore(const PairedReadsDatastore & datastore) {paired_datastores.emplace_back(&datastore);}
    void add_datastore(const LongReadsRecruiter & datastore) {long_recruiters.emplace_back(&datastore);}
    //walk out from a node
    SequenceDistanceGraphPath walk_out(sgNodeID_t n);

    SequenceDistanceGraphPath walk_out_in_order(sgNodeID_t n, bool use_pair=true, bool collapse_pair=true, bool verbose=false);
    WorkSpace & ws;
    std::vector<const PairedReadsDatastore *> paired_datastores;
    std::vector<const LongReadsRecruiter *> long_recruiters;
    static const std::string logo;
};

