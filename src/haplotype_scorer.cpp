//
// Created by Katie Barr (EI) on 14/11/2017.
//

#include "haplotype_scorer.h"
double avg(std::vector<int> v){
    if (v.size() > 0) {
        return std::accumulate(v.begin(), v.end(), 0LL) / v.size();
    }
    return 0.0;
}

double stdev(std::vector<int> v, double mean){
    if (v.size() > 0) {
        double res = 0;
        for (auto i: v) {
            res += std::pow(i - mean, 2);
        }
        return std::pow(res / v.size(), 0.5);
    }
    return 0.0;
}



void print_vector(std::vector<std::string> vec){
    for (auto a: vec){
        std::cout << a << " ";
    }
    std::cout << std::endl;
}

void print_vector_ids(std::vector<sgNodeID_t > vec, std::map<sgNodeID_t , std::string > id_to_contig_name){
    for (auto a: vec){
        std::cout << id_to_contig_name[a] << " ";
    }
    std::cout << std::endl;
}

void print_int_vector(std::vector<int> vec){
    for (auto a: vec){
        std::cout << a << " ";
    }
    std::cout << std::endl;
}

void print_pair_int_vector(std::vector<std::pair<int, int> > vec){
    for (auto a: vec){
        std::cout << std::get<0>(a) << " " << std::get<1>(a);
    }
    std::cout << std::endl;
}
HaplotypeScorer::HaplotypeScorer(SequenceGraph & sg, std::string fasta_filename): sg(sg), mapper(PairedReadMapper(sg)){

}

void HaplotypeScorer::decide_barcode_haplotype_support(){

    int support;
    int haplotypes_supported = 0;

    std::cout << "Calculating barcode haplotype" <<std::endl;
    for (auto &mapping:barcode_node_mappings) {
        if (mapping.second.size() > 1) {
            std::vector<sgNodeID_t> nodes;
            std::vector<int> scores;
            for (auto n: mapping.second) {
                nodes.push_back(n.first);
                scores.push_back(n.second);
            }
            if (*std::max_element(scores.begin(), scores.end()) > 1) {
                //  for i, haplotype in enumerate(self.list_of_possible_haplotypes)
                for (int i = 0; i < haplotype_ids.size(); i++) {
                    std::vector<sgNodeID_t> nodes_in_haplotype;
                    std::vector<sgNodeID_t> h;
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
        // try alternative way just looping over each node mapped to and each hap for that node
        if (mapping.second.size() > 1) {
            for (auto n:mapping.second) {
                // won't be idemtical to before because allows scores to all be 1
                // see if this agrees well with other
                for (auto h: node_id_haplotype_index_map[n.first]) {
                    barcode_haplotype_mappings2[mapping.first][h] += n.second;

                }
            }
        }

    }
    bool no_matched = true;
    for (auto m:barcode_haplotype_mappings){

        for (auto n:m.second){
            if (barcode_haplotype_mappings2[m.first][n.second] != barcode_haplotype_mappings[m.first][n.second]){
                no_matched = false;
                break;
            }
        }
    }
    if (!no_matched){
        std::cout << "twu methods differ " << std::endl;
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
                id_to_contig_name[node] = n;
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
    std::cout << "Mapped " << mapper.read_to_node.size() << " reads to " <<  mapper.reads_in_node.size() << "nodes" << std::endl;
    std::cout << "NOw counting barcode votes... " << std::endl;
    int counter = 0;
    for (auto &node:haplotype_nodes){
        //std::cout << "Node: " <<node <<std::endl;
        for (auto &mapping:mapper.reads_in_node[node>0?node:-node]){
            std::string barcode = mapper.read_ids_to_tags[mapping.read_id];
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
    std::vector<int> haplotype_support_vals;
    std::vector<int> haplotype_not_support_vals;
    std::vector<int> haplotype_overall_support_vals;
    for (int i = 0; i < haplotype_ids.size(); i++) {
        haplotype_support_vals.push_back(haplotype_support[i]);
        haplotype_not_support_vals.push_back(haplotype_not_support[i]);
        haplotype_overall_support_vals.push_back(haplotype_overall_support[i]);

    }
    std::vector<std::pair<int, int> > supports;
    std::vector<std::pair<int, int> > overall_supports;
    std::vector<std::pair<std::pair<int, int>, int > > pair_supports;
    std::vector<std::pair<std::pair<int, int>, int > > pair_overall_supports;
    for (int i=0; i< haplotype_ids.size(); i++){
        supports.push_back(std::make_pair(i, haplotype_support[i]));
        overall_supports.push_back(std::make_pair(i, haplotype_support[i]));
    }
    std::sort(supports.begin(), supports.end(), [](auto &left, auto &right) {
        return left.second < right.second;
    });
    std::sort(overall_supports.begin(), overall_supports.end(), [](auto &left, auto &right) {
        return left.second < right.second;
    });
    // when not rushing, get rid of al this repetition
    auto support_max = std::max_element(haplotype_support_vals.begin(), haplotype_support_vals.end());
    auto support_mean = avg(haplotype_support_vals);
    auto overall_support_max = std::max_element(haplotype_overall_support_vals.begin(),
                                                haplotype_overall_support_vals.end());
    auto overall_support_mean = avg(haplotype_support_vals);

    auto support_stdev = stdev(haplotype_support_vals, support_mean);

    auto overall_stdev = stdev(haplotype_overall_support_vals, overall_support_mean);

    //print_summary(outfile, supports, overall_supports, haplotype_support_vals, haplotype_overall_support_vals, haplotype_not_support_vals);
    if (hap_pair_support.size() > 0 && hap_pair_support_total_score.size() > 0) {
        std::vector<int> hap_pair_not_support_values;
        std::vector<int> hap_pair_support_values;
        for (auto h : hap_pair_support) {
            hap_pair_support_values.push_back(h.second);
            pair_supports.push_back(std::make_pair(h.first, h.second));
        }
        for (auto h : hap_pair_not_support) {
            hap_pair_not_support_values.push_back(h.second);
        }
        std::vector<int> hap_pair_support_total_score_values;
        for (auto h : hap_pair_support_total_score) {
            hap_pair_support_total_score_values.push_back(h.second);
            pair_overall_supports.push_back(std::make_pair(h.first, h.second));
        }
        std::sort(pair_supports.begin(), pair_supports.end(), [](auto &left, auto &right) {
            return left.second < right.second;
        });
        std::sort(pair_overall_supports.begin(), pair_overall_supports.end(), [](auto &left, auto &right) {
            return left.second < right.second;
        });

        //print_pair_summary(outfile, pair_supports, pair_overall_supports, hap_pair_support_values, hap_pair_support_total_score_values, hap_pair_not_support_values);

        // get winners
        std::vector<int> support_winner;
        std::vector<int> overall_support_winner;
        for (int h = 0; h < haplotype_ids.size(); h++) {
            if (haplotype_support[h] == *support_max) {
                support_winner.push_back(h);
            }
            if (haplotype_overall_support[h] == *overall_support_max) {
                overall_support_winner.push_back(h);
            }
        }
        std::vector<std::pair<int, int> > pair_support_winner;
        std::vector<std::pair<int, int> > pair_overall_support_winner;
        auto overall_pair_support_max = std::max_element(hap_pair_support_total_score_values.begin(),
                                                         hap_pair_support_total_score_values.end());

        auto pair_support_max = std::max_element(hap_pair_support_values.begin(), hap_pair_support_values.end());

        for (auto h: hap_pair_support) {
            if (h.second == *pair_support_max) {
                pair_support_winner.push_back(h.first);
            }
        }
        for (auto h: hap_pair_support_total_score) {
            if (h.second == *overall_pair_support_max) {
                pair_overall_support_winner.push_back(h.first);

            }
        }
        std::cout << "Support winner: ";
        print_int_vector(support_winner);
        std::cout << "overall SUpport winner: ";
        print_int_vector(overall_support_winner);
        std::cout << "pair SUpport winner: ";
        print_pair_int_vector(pair_support_winner);
        std::cout << "pair overall SUpport winner: ";
        print_pair_int_vector(pair_overall_support_winner);
        print_vector_ids(haplotype_ids[support_winner[0]], id_to_contig_name);
        int winner = *std::max_element(haplotype_support, haplotype_support + number_haplotypes);
        std::cout << "Winner: " << winner << std::endl;
        for (auto h: haplotype_ids[winner]) {
            std::cout << id_to_contig_name[h] << " ";
        }
        std::cout << std::endl;
    }
    return 0;
}


int HaplotypeScorer::score_haplotypes2() {
    auto number_haplotypes = haplotype_ids.size();
    std::cout << "Finding most supported of " << number_haplotypes<< " possible haplotypes"<<std::endl;
    //initialize score arrays- index is haplotype index
    int haplotype_support[number_haplotypes] = {0};
    int haplotype_not_support[number_haplotypes] = {0};
    int haplotype_overall_support[number_haplotypes] = {0};
    std::map<std::pair<int, int>, int> hap_pair_not_support;
    std::map<std::pair<int, int>, int> hap_pair_support;
    std::map<std::pair<int, int>, int> hap_pair_support_total_score;

    std::map<int, std::map<std::string, int > > haplotype_barcode_agree2;
    std::map<int, std::map<std::string, int > > haplotype_barcode_disagree2;
    std::string barcode;
    for (auto &bm: barcode_haplotype_mappings2) {
        barcode = bm.first;
        std::vector<int> winners = winner_for_barcode(barcode); // ideally should be length 1
        for (auto winner:winners) {
            int pair = haplotype_ids.size() - 1 - winner;
            haplotype_support[winner] += 1;
            hap_pair_support[std::make_pair(winner, pair)] += 1;
            haplotype_barcode_agree2[winner][barcode] += bm.second[winner];
            haplotype_barcode_disagree2[winner][barcode] += bm.second[pair];
        }
        for (int hap = 0; hap < number_haplotypes / 2; hap++) {
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
    std::vector<int> haplotype_support_vals;
    std::vector<int> haplotype_not_support_vals;
    std::vector<int> haplotype_overall_support_vals;
    for (int i = 0; i < haplotype_ids.size(); i++) {
        haplotype_support_vals.push_back(haplotype_support[i]);
        haplotype_not_support_vals.push_back(haplotype_not_support[i]);
        haplotype_overall_support_vals.push_back(haplotype_overall_support[i]);

    }
    std::vector<std::pair<int, int> > supports;
    std::vector<std::pair<int, int> > overall_supports;
    std::vector<std::pair<std::pair<int, int>, int > > pair_supports;
    std::vector<std::pair<std::pair<int, int>, int > > pair_overall_supports;
    for (int i=0; i< haplotype_ids.size(); i++){
        supports.push_back(std::make_pair(i, haplotype_support[i]));
        overall_supports.push_back(std::make_pair(i, haplotype_support[i]));
    }
    std::sort(supports.begin(), supports.end(), [](auto &left, auto &right) {
        return left.second < right.second;
    });
    std::sort(overall_supports.begin(), overall_supports.end(), [](auto &left, auto &right) {
        return left.second < right.second;
    });
    // when not rushing, get rid of al this repetition
    auto support_max = std::max_element(haplotype_support_vals.begin(), haplotype_support_vals.end());
    auto support_mean = avg(haplotype_support_vals);
    auto overall_support_max = std::max_element(haplotype_overall_support_vals.begin(),
                                                haplotype_overall_support_vals.end());
    auto overall_support_mean = avg(haplotype_support_vals);

    auto support_stdev = stdev(haplotype_support_vals, support_mean);

    auto overall_stdev = stdev(haplotype_overall_support_vals, overall_support_mean);

    //print_summary(outfile, supports, overall_supports, haplotype_support_vals, haplotype_overall_support_vals, haplotype_not_support_vals);
    if (hap_pair_support.size() > 0 && hap_pair_support_total_score.size() > 0) {
        std::vector<int> hap_pair_not_support_values;
        std::vector<int> hap_pair_support_values;
        for (auto h : hap_pair_support) {
            hap_pair_support_values.push_back(h.second);
            pair_supports.push_back(std::make_pair(h.first, h.second));
        }
        for (auto h : hap_pair_not_support) {
            hap_pair_not_support_values.push_back(h.second);
        }
        std::vector<int> hap_pair_support_total_score_values;
        for (auto h : hap_pair_support_total_score) {
            hap_pair_support_total_score_values.push_back(h.second);
            pair_overall_supports.push_back(std::make_pair(h.first, h.second));
        }
        std::sort(pair_supports.begin(), pair_supports.end(), [](auto &left, auto &right) {
            return left.second < right.second;
        });
        std::sort(pair_overall_supports.begin(), pair_overall_supports.end(), [](auto &left, auto &right) {
            return left.second < right.second;
        });

        //print_pair_summary(outfile, pair_supports, pair_overall_supports, hap_pair_support_values, hap_pair_support_total_score_values, hap_pair_not_support_values);

        // get winners
        std::vector<int> support_winner;
        std::vector<int> overall_support_winner;
        for (int h = 0; h < haplotype_ids.size(); h++) {
            if (haplotype_support[h] == *support_max) {
                support_winner.push_back(h);
            }
            if (haplotype_overall_support[h] == *overall_support_max) {
                overall_support_winner.push_back(h);
            }
        }
        std::vector<std::pair<int, int> > pair_support_winner;
        std::vector<std::pair<int, int> > pair_overall_support_winner;
        auto overall_pair_support_max = std::max_element(hap_pair_support_total_score_values.begin(),
                                                         hap_pair_support_total_score_values.end());

        auto pair_support_max = std::max_element(hap_pair_support_values.begin(), hap_pair_support_values.end());

        for (auto h: hap_pair_support) {
            if (h.second == *pair_support_max) {
                pair_support_winner.push_back(h.first);
            }
        }
        for (auto h: hap_pair_support_total_score) {
            if (h.second == *overall_pair_support_max) {
                pair_overall_support_winner.push_back(h.first);

            }
        }
        std::cout << "Support winner: ";
        print_int_vector(support_winner);
        std::cout << "overall SUpport winner: ";
        print_int_vector(overall_support_winner);
        std::cout << "pair SUpport winner: ";
        print_pair_int_vector(pair_support_winner);
        std::cout << "pair overall SUpport winner: ";
        print_pair_int_vector(pair_overall_support_winner);
        print_vector_ids(haplotype_ids[support_winner[0]], id_to_contig_name);
        // this doesn't work cos hap dupport is a ddictionary
        int winner = *std::max_element(haplotype_support, haplotype_support + number_haplotypes);
        std::cout << "Winner method 2: " << winner << std::endl;
        for (auto h: haplotype_ids[winner]) {
            std::cout << id_to_contig_name[h] << " ";
        }
        std::cout << std::endl;

    }
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