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
    int counter = 0;
    std::string res;
    // no obvious way to access graph nodes directly, so store bubble names;
    std::vector<std::vector<std::string> > bubbles;
    while (std::getline(infile, line)) {
        std::istringstream(line) >> res;
        if (res == "next" or counter == max_degree){
            std::cout << "next" << std::endl;
            std::vector<std::string> next_bubble;
            for (int j=0; j < counter; j++){
                std::cout << "addng: " << fields[j] << std::endl;
                next_bubble.push_back(fields[j]);
            }

            counter = 0;
            bubbles.push_back(next_bubble);
            for (auto b: bubbles.back()){
                std::cout << "b: " << b << std::endl;
            }

            continue;
        } else {
            std::cout << "res: " << res << std::endl;
            fields[counter] = res;
            counter += 1;

        }

    }
    std::cout << "Bubble \n";
    for (auto b:bubbles){
        for (auto c:b){
            std::cout << c << " ";
        }
        std::cout << "\n";
    }
    std::cout << "\n";
};
void HaplotypeScorer::count_barcode_votes(){};
void HaplotypeScorer::score_haplotypes(){};