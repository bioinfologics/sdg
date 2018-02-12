//
// Created by Katie Barr (EI) on 23/11/2017.
//

#ifndef SG_PHASESCAFFOLDER_H
#define SG_PHASESCAFFOLDER_H

#include <string>


#include "sglib/SequenceGraph.h"
#include <sglib/PairedReadMapper.h>
#include <sglib/HaplotypeScorer.hpp>
/**
 *
 *
 * @brief
 * Find haplotypes for each component and join into phased blocks
 *
 */

class PhaseScaffolder {
public:
    PhaseScaffolder(SequenceGraph &);
    SequenceGraph & sg;

    void output_bubbles(std::string);
    void phase_components();
    void load_mappings(std::string , std::string , std::string, uint64_t );
    PairedReadMapper mapper;

};


#endif //SG_PHASESCAFFOLDER_H
