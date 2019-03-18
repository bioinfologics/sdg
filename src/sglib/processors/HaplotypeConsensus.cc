//
// Created by Bernardo Clavijo (EI) on 2019-02-27.
//

#include "HaplotypeConsensus.hpp"

std::string HaplotypeConsensus::run_consensus() {
    add_non_specific_nodes();
    return consensus_sequence();
}


std::string HaplotypeConsensus::consensus_sequence() {
    return "";
}

void HaplotypeConsensus::use_long_reads_from_file(std::string filename) {

}