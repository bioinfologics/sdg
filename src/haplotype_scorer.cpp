//
// Created by Katie Barr (EI) on 14/11/2017.
//

#include "haplotype_scorer.h"

HaplotypeScorer::HaplotypeScorer(SequenceGraph & sg, std::string fasta_filename): sg(sg), mapper(PairedReadMapper(sg)){

}

void HaplotypeScorer::decide_barcode_haplotype_support(){

    int support;
    int haplotypes_supported = 0;
    std::cout << "Calculating barcode haplotype support for" <<std::endl;
    for (auto &mapping:barcode_node_mappings){
        if (mapping.second.size() > 1){
            std::vector<sgNodeID_t > nodes;
            std::vector<int> scores;
            for (auto n: mapping.second){
                nodes.push_back(n.first);
                scores.push_back(n.second);
            }
            if (*std::max_element(scores.begin(), scores.end())> 1) {
                //  for i, haplotype in enumerate(self.list_of_possible_haplotypes)
                for (int i = 0; i < haplotype_ids.size(); i++) {
                    std::vector<sgNodeID_t > nodes_in_haplotype;
                    std::vector<sgNodeID_t > h;
                    h = haplotype_ids[i];
                    // find all nods in each haplotype that this barcode maps to
                    for (auto n1: nodes) {
                        //edges_in_haplotype = [e for e in haplotype if e in edges]
                        if (std::find(h.begin(), h.end(), n1) != h.end()) {
                            nodes_in_haplotype.push_back(n1);
                        }
                    }
                    // somewhat arbitrary rule to decide if the barcode supports a haplotype enough
                    // if len(edges_in_haplotype)>= len(edges)/2 and len(edges_in_haplotype) > 1:
                    if (nodes_in_haplotype.size() >= (nodes.size() / 2) && nodes_in_haplotype.size() > 1) {
                        support = 0;
                        for (auto a: nodes_in_haplotype) {
                            support += mapping.second[a];
                        }
                        barcode_haplotype_mappings[mapping.first][i] = support;
                        support = 0;
                        haplotypes_supported += 1;
                    } else {
                        unused_barcodes.push_back(mapping.first);
                    }
                }
            }

        } else {
            unused_barcodes.push_back(mapping.first);
        }
        //std::cout << "barcode " << mapping.first << " supports " << haplotypes_supported << std::endl;

        haplotypes_supported = 0;
        // try alternative way just looping over each edge mapped to and each hap in that edge
        /*for (auto n:mapping.second){

        }*/
    }
    std::cout << "Calculated haplotype support for each barcode, " << barcode_haplotype_mappings.size() <<  std::endl;

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
    std::vector<std::string> fields;
    std::cout << "Loading haplotypes " << haplotype_filename << std::endl;
    std::string res;
    // no obvious way to access graph nodes directly, so store  names;
    std::vector<std::vector<sgNodeID_t> > haplotypes;
    while (std::getline(infile, line)) {
        std::istringstream(line) >> res;
        if (res == "next" ){
            std::vector<sgNodeID_t> next_haplotype;
            for (auto n:fields){
                sgNodeID_t node = sg.oldnames_to_ids[n];
                next_haplotype.push_back(node);
                haplotype_nodes.insert(node);
                node_id_haplotype_index_map[node].push_back(haplotypes.size());
                fields.clear();
            }

            haplotypes.push_back(next_haplotype);

            continue;
        } else {
            fields.push_back(res);

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
    mapper.map_reads(r1_filename, r2_filename, prm10x);
    // loop over all mappings to get reads mapped to haplotype nodes - can only get kmers from read in node,
    /*std::vector<int> reads_mapped_to_het_nodes;
    for (int read_index= 0; read_index < mapper.read_to_node.size(); read_index++) {
        std::string barcode = mapper.read_to_tag[read_index];
        sgNodeID_t node = mapper.read_to_node[read_index];
        // if read mapped succesfully
         if (node > 0){
        //    barcode_node_mappings[barcode][node]  += mapper.
             if (std::find(haplotype_nodes.begin(), haplotype_nodes.end(), node) != haplotype_nodes.end()){
                 reads_mapped_to_het_nodes.push_back(read_index);

             }
        }

    }*/
    std::cout << "Mapped " << mapper.read_to_node.size() << " reads to " <<  mapper.reads_in_node.size() << "nodes" << std::endl;
    std::cout << "NOw counting barcode votes... " << std::endl;
    int counter = 0;
    for (auto &node:haplotype_nodes){
        std::cout << "Node: " <<node <<std::endl;
        for (auto &mapping:mapper.reads_in_node[node>0?node:-node]){
            std::string barcode = mapper.read_ids_to_tags[mapping.read_id];
            std::cout << "Barcode: " << barcode << std::endl;
             barcode_node_mappings[barcode][node] += mapping.unique_matches;
            counter += 1;
        }
    }
    std::cout << "counted " << counter << " votes for " << barcode_node_mappings.size() << "barcodes";
};

int HaplotypeScorer::score_haplotypes() {
    auto number_haplotypes = haplotype_ids.size();
    std::cout << "Finding most supported of " << number_haplotypes<< " possible haplotypes"<<std::endl;
    //initialize score arrays- index is haplotype index
    int haplotype_support[number_haplotypes] = {0};
    int haplotype_not_support[number_haplotypes] = {0};
    int haplotype_overall_support[number_haplotypes] = {0};
    std::map<std::pair<int, int>, int> hap_pair_not_support;
    std::map<std::pair<int, int>, int> hap_pair_support;
    std::map<std::pair<int, int>, int> hap_pair_support_total_score;
    std::string barcode;
    for (auto &bm: barcode_haplotype_mappings) {
        barcode = bm.first;
        std::vector<int> winners = winner_for_barcode(barcode); // ideally should be length 1
        for (auto winner:winners){
            int pair = haplotype_ids.size() - 1 - winner;
            haplotype_support[winner] += 1;
            hap_pair_support[std::make_pair(winner, pair)] += 1;
            haplotype_barcode_agree[winner][barcode] += bm.second[winner];
            haplotype_barcode_disagree[winner][barcode] += bm.second[pair];
        }
        for (int hap = 0; hap < number_haplotypes/ 2; hap++) {
            // pair = len(self.list_of_possible_haplotypes) -1 -haplotype
            auto pair = number_haplotypes - 1 - hap;
            if (bm.second.find(hap) != bm.second.end()) {
                haplotype_overall_support[hap] += bm.second[hap];
                hap_pair_support_total_score[std::make_pair(hap, pair)] += bm.second[hap];
            }

            if (bm.second.find(pair) != bm.second.end()) {
                haplotype_overall_support[pair] += bm.second[pair];
                hap_pair_support_total_score[std::make_pair(hap, pair)] += bm.second[pair];
            }
            if (bm.second.find(hap) == bm.second.end()) {
                haplotype_not_support[hap] += 1;
            }
            if (bm.second.find(pair) == bm.second.end()) {
                haplotype_not_support[pair] += 1;
            }
            if (bm.second.find(hap) == bm.second.end() and bm.second.find(pair) == bm.second.end()) {
                hap_pair_not_support[std::make_pair(hap, pair)] += 1;

            }
        }
    }
    int winner = *std::max_element(haplotype_support, haplotype_support+number_haplotypes);
    std::cout << "Winner: " << winner << std::endl;
    for (auto h: haplotype_ids[winner]){
        std::cout << h << " ";
    }
    std::cout << std::endl;
    return 0;
}


std::vector<int>  HaplotypeScorer::winner_for_barcode(std::string barcode){
    int max=0;
    std::vector<int> winners;
    for (auto h:barcode_haplotype_mappings[barcode]){
        if (h.second > max){
            max = h.second;
        }
    }
    //TODO: DECIDE CRITERIA FOR MINIMUM SUPPORT
    for (auto h:barcode_haplotype_mappings[barcode]){
        if (h.second == max){
            winners.push_back(h.first);
        }
    }
    return winners;
}