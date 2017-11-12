//
// Created by Bernardo Clavijo (EI) on 10/11/2017.
//

#ifndef SG_SCAFFOLDER_HPP
#define SG_SCAFFOLDER_HPP

#include "SequenceGraph.hpp"
#include "PairedReadMapper.hpp"

class Scaffolder {

public:
    Scaffolder(SequenceGraph &_sg, std::vector<PairedReadMapper> & _rms) : sg(_sg),rmappers(_rms){};

    //1: find and expand canonical repeats
    void find_canonical_repeats();

    uint64_t count_reads_linking(sgNodeID_t source, sgNodeID_t dest);

    //2: find unsatisfied connections
    //find read support breakpoints
    //3: path finder?


    SequenceGraph & sg;
    std::vector<PairedReadMapper> &rmappers;




};


#endif //SG_SCAFFOLDER_HPP
