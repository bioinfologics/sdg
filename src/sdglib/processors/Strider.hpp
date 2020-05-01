//
// Created by Bernardo Clavijo (EI) on 01/05/2020.
//

#pragma once

#include <sdglib/workspace/WorkSpace.hpp>

class Strider {
public:

    Strider(WorkSpace & _ws):ws(_ws){};
    SequenceDistanceGraphPath walk_out(sgNodeID_t n);
    WorkSpace & ws;
};

