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
#include <numeric>

#include "sglib/SequenceGraph.hpp"
#include <sglib/PairedReadMapper.hpp>

class HaplotypeScorer{

public:
     HaplotypeScorer(SequenceGraph &);
    // functions we will need:
    void find_possible_haplotypes(int, std::string);
    void load_haplotypes(std::string, int);

    void count_barcode_votes(std::string, std::string);
    int score_haplotypes();
    int score_haplotypes2();

    std::map<std::string, std::map<sgNodeID_t , int> > barcode_node_mappings;
    void decide_barcode_haplotype_support();
    std::map<std::string, std::map< int, int > > barcode_haplotype_mappings;
    std::map<std::string, std::map< int, int > > barcode_haplotype_mappings2;

private:
    SequenceGraph & sg;
    PairedReadMapper mapper;
    // each het node
    std::set<sgNodeID_t > haplotype_nodes;
    // each possible hsplotype
    std::vector<std::vector<sgNodeID_t> > haplotype_ids;

    std::vector <std::string> unused_barcodes;
    std::map<sgNodeID_t , std::vector<int> > node_id_haplotype_index_map;
    std::map<sgNodeID_t , std::string > id_to_contig_name;

            std::vector<int>  winner_for_barcode(std::string barcode);

    std::map<int, std::map<std::string, int > > haplotype_barcode_agree;
    std::map<int, std::map<std::string, int > > haplotype_barcode_disagree;

};
#endif //SG_HAPLOTYPE_SCORER_H
