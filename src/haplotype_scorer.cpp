//
// Created by Katie Barr (EI) on 14/11/2017.
//

#include "haplotype_scorer.h"

HaplotypeScorer::HaplotypeScorer(SequenceGraph sg, std::string fasta_filename): sg(sg), mapper(sg){

}

void HaplotypeScorer::find_possible_haplotypes(int max_degree, std::string bubble_contig_filename){
    std::ifstream infile(bubble_contig_filename);
    std::string line;
    std::string fields[max_degree];
    std::cout << "Loading bubbles " << bubble_contig_filename << std::endl;

    while (std::getline(infile, line)) {
        std::string previous = "";
        std::string current;
        std::string current2;
        std::cout << line << std::endl;
        for (int i=0; i < max_degree; i++){
            std::istringstream(line) >> current;
            std::istringstream(line) >> current2;
            std::cout << "i: " << i << "Current: " << current  << " Current2: " << current2 << " previous " << previous << std::endl;
            if (current!= previous) {
                fields[i] = current;
                previous = fields[i];
            } else {
                //  get same contig name twice when its run past last column
                break;
            }
        }
        std::cout << "Bubble \n";
        for (int i=0; i < max_degree; i++){
            std::cout << fields[i] << " ";
        }
        std::cout << "\n";

    }
};
void HaplotypeScorer::count_barcode_votes(){};
void HaplotypeScorer::score_haplotypes(){};