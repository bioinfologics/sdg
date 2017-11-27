//
// Created by Katie Barr (EI) on 23/11/2017.
//

#include "PhaseScaffolder.h"



PhaseScaffolder::PhaseScaffolder(SequenceGraph & sg): sg(sg), mapper(sg){
}

void PhaseScaffolder::load_mappings(std::string r1_filename, std::string r2_filename, std::string fasta_filename, uint64_t max_mem_gb){

    mapper.map_reads(r1_filename, r2_filename, fasta_filename, prm10x, max_mem_gb);
    std::cout << "Mapped " << mapper.read_to_node.size() << " reads to " <<  mapper.reads_in_node.size() << "nodes" << std::endl;
}


void PhaseScaffolder::output_bubbles(std::string bubble_filename) {

        std::ofstream out(bubble_filename);
//find each component of gfa
    auto components = sg.connected_components();
// this finds 2 components for test graph...
    std::cout << "Found " << components.size() << " connected components " << std::endl;
    int counter = 0;
    for (auto component:components) {
        auto bubbles = sg.find_bubbles(component);
        if(bubbles.size() > 1){
            for (auto bubble:bubbles)
            {
                for (auto bubble_c:bubble){
                    out << ">"<<sg.oldnames[bubble_c] << "_" << counter << std::endl << sg.nodes[bubble_c].sequence << std::endl;
                }

            }

        }
        counter +=1;

    }
}

void PhaseScaffolder::phase_components() {

//find and phase each component of gfa
    auto components = sg.connected_components();
// this finds 2 components for test graph...
    std::cout << "Found " << components.size() << " connected components " << std::endl;
    for (auto component:components) {
        HaplotypeScorer hs;

// should
        auto bubbles = sg.find_bubbles(component);
        hs.find_possible_haplotypes(bubbles);
        std::cout << "mapper.reads_in_node.size()  " << mapper.reads_in_node.size() <<std::endl;

        hs.count_barcode_votes(mapper);
        hs.decide_barcode_haplotype_support();
// now have mappings and barcode support
        if (hs.barcode_haplotype_mappings.size() > 0) {
            hs.score_haplotypes(sg.oldnames);
            std::cout << "scored haplotypes " << hs.success << std::endl;
        }
    }
}