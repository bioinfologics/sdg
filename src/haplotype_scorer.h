//
// Created by Katie Barr (EI) on 14/11/2017.
//

#ifndef SG_HAPLOTYPE_SCORER_H
#define SG_HAPLOTYPE_SCORER_H

#include "sglib/SequenceGraph.hpp"
#include <sglib/PairedReadMapper.hpp>

class HaplotypeScorer{

public:
     HaplotypeScorer(SequenceGraph, std::string);

private:
    SequenceGraph sg;
    std::string mapping_filename;
    PairedReadMapper mapper;
};
#endif //SG_HAPLOTYPE_SCORER_H
