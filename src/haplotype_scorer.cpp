//
// Created by Katie Barr (EI) on 14/11/2017.
//

#include "haplotype_scorer.h"

HaplotypeScorer::HaplotypeScorer(SequenceGraph sg, std::string fasta_filename): sg(sg), mapper(sg){

}

std::vector<std::vector<std::string> > HaplotypeScorer::load_bubble_file(std::string bubble_contig_filename, int max_degree){
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
            std::vector<std::string> next_bubble;
            for (int j=0; j < counter; j++){
                next_bubble.push_back(fields[j]);
            }

            counter = 0;
            bubbles.push_back(next_bubble);

            continue;
        } else {
            fields[counter] = res;
            counter += 1;

        }

    }
    std::cout << "Loaded " << bubbles.size() << "bubbles \n";
    return  bubbles;
}

void HaplotypeScorer::find_possible_haplotypes(int max_degree, std::string bubble_contig_filename){

    auto bubbles = load_bubble_file(bubble_contig_filename, max_degree);
};

void  HaplotypeScorer::load_haplotypes(std::string haplotype_filename, int degree = 2){
    std::ifstream infile(haplotype_filename);
    std::string line;
    std::string fields[degree];
    std::cout << "Loading haplotypes " << haplotype_filename << std::endl;
    int counter = 0;
    std::string res;
    // no obvious way to access graph nodes directly, so store  names;
    std::vector<std::vector<sgNodeID_t> > haplotypes;
    while (std::getline(infile, line)) {
        std::istringstream(line) >> res;
        if (res == "next" or counter == degree){
            std::vector<sgNodeID_t> next_haplotype;
            for (int j=0; j < counter; j++){
                next_haplotype.push_back(sg.oldnames_to_ids[fields[j]]);
            }

            counter = 0;
            haplotypes.push_back(next_haplotype);

            continue;
        } else {
            fields[counter] = res;
            counter += 1;

        }
        //TODO: sanity check, no ids should be 0, all haps should be same length


    }
    haplotype_ids =  haplotypes;
    for (auto h: haplotypes){
        for (auto c: h) {
            std::cout << c << " ";

        }
        std::cout << std::endl << "next:";
    }
    std::cout << "Loaded " << haplotypes.size() << "haplotypes \n";
}

void HaplotypeScorer::count_barcode_votes(std::string r1_filename, std::string r2_filename){
    mapper.map_reads(r1_filename, r2_filename);
};
void HaplotypeScorer::score_haplotypes(){};