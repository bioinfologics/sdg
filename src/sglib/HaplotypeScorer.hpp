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

#include "sglib/SequenceGraph.h"
#include <sglib/PairedReadMapper.h>

class HaplotypeScorer{

public:
    // functions we will need:
    void find_possible_haplotypes(std::vector<std::vector<sgNodeID_t >>);
    void load_haplotypes(std::string, int);

    void count_barcode_votes(PairedReadMapper &);
    int score_haplotypes(std::vector<std::string> );

    std::map<prm10xTag_t, std::map<sgNodeID_t , int> > barcode_node_mappings;
    void decide_barcode_haplotype_support();
    std::map<prm10xTag_t, std::map< int, int > > barcode_haplotype_mappings;
    std::map<prm10xTag_t, std::map< int, int > > barcode_haplotype_mappings2;
    bool success = false; // if doing partial success replace with enum

private:

    // each het node
    std::set<sgNodeID_t > haplotype_nodes;
    // each possible hsplotype
    std::vector<std::vector<sgNodeID_t> > haplotype_ids;

    std::vector <prm10xTag_t> unused_barcodes;
    std::map<sgNodeID_t , std::vector<int> > node_id_haplotype_index_map;
    std::map<sgNodeID_t , std::string > id_to_contig_name;

    std::vector<int>  winner_for_barcode(prm10xTag_t barcode);

    std::map<int, std::map<prm10xTag_t, int > > haplotype_barcode_agree;
    std::map<int, std::map<prm10xTag_t, int > > haplotype_barcode_disagree;

};
#endif //SG_HAPLOTYPE_SCORER_H
