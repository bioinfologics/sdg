//
// Created by Katie Barr (EI) on 14/11/2017.
//

#ifndef SG_HAPLOTYPE_SCORER_H
#define SG_HAPLOTYPE_SCORER_H
#include <sstream>
#include <iostream>
#include <fstream>
#include <sstream>
#include <istream>
#include <string>
#include "sglib/SequenceGraph.hpp"
#include <sglib/PairedReadMapper.hpp>

class HaplotypeScorer{

public:
     HaplotypeScorer(SequenceGraph, std::string);
    // functions we will need:
    void find_possible_haplotypes(int, std::string);
    void load_haplotypes(std::string, int);

    void count_barcode_votes(std::string, std::string);
    void score_haplotypes();
    std::map<std::string, std::map<sgNodeID_t , int> > barcode_node_mappings;

private:
    SequenceGraph sg;
    std::string fasta_filename;
    PairedReadMapper mapper;
    // each het node
    std::set<sgNodeID_t > haplotype_nodes;
    std::vector<std::vector<std::string> > load_bubble_file(std::string , int );
    // each possible hsplotype
    std::vector<std::vector<sgNodeID_t> > haplotype_ids;

};
#endif //SG_HAPLOTYPE_SCORER_H
