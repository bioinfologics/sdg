//
// Created by Bernardo Clavijo (EI) on 2019-02-27.
//

#include "HaplotypeConsensus.hpp"
#include <sglib/utilities/most_common_helper.hpp>

std::string HaplotypeConsensus::consensus_sequence() {
    std::string consensus;
    std::vector<sgNodeID_t > current_path;

    for (uint32_t ni = 0; ni < backbone_filled_path.size(); ++ni) {
        const auto &n = backbone_filled_path[ni];
        if (n!=0) {
            current_path.emplace_back(n);
        }
        if (n == 0 or ni==backbone_filled_path.size()-1) {
            consensus+=SequenceGraphPath(ldg.sg, current_path).get_sequence(); // TODO: CAREFUL!!! this is marked as failing!!
            current_path.clear();
        }
        if (n == 0 and ni < backbone_filled_path.size()-1){
            int dist=200;
            for (const auto &fwl: ldg.get_fw_links(backbone_filled_path[ni-1])) {
                if (fwl.dest == backbone_filled_path[ni+1]) dist=fwl.dist;
            }
            consensus += std::string(std::max(dist,100), 'N');
        }
    }
    return consensus;
}

void HaplotypeConsensus::use_long_reads_from_file(std::string filename) {
    std::ifstream reads(filename);
    if (!reads.good()) {
        std::cerr << "Error opening: " << filename << ".\n" << std::strerror(errno) << std::endl;
        return;
    }
    std::string line, id, DNA_sequence;
    uint32_t rid(0);
    while (std::getline(reads, line).good()) {
        if (line[0] == '>') {
            id = line.substr(1);
            if (!DNA_sequence.empty()) {
                read_seqs[rid]=DNA_sequence;
            }
            rid = std::stoul(id);
            DNA_sequence.clear();
        } else if (line[0] != '>') {
            DNA_sequence += line;
        }
    }

    return;
}

void HaplotypeConsensus::orient_read_path(uint32_t rid) {
    std::set<sgNodeID_t > l(backbone.cbegin(), backbone.cend());

    const auto &forward_path = ws.long_read_mappers[0].read_paths[rid];
    std::cout << "Original read path: " << std::endl;
    std::copy(forward_path.cbegin(), forward_path.cend(), std::ostream_iterator<sgNodeID_t>(std::cout, ", "));
    std::cout << std::endl;

    std::vector<sgNodeID_t> reversed_path;
    uint32_t count_reversed(0);
    uint32_t count_forward(0);
    // Look at the path reversed and with opposite signs
    std::transform(forward_path.crbegin(), forward_path.crend(), std::back_inserter(reversed_path), [](const sgNodeID_t &n) {return -n;});
    std::for_each(reversed_path.cbegin(), reversed_path.cend(),
            [&l,&count_reversed](const sgNodeID_t &n){l.find(n) != l.cend() ? ++count_reversed:count_reversed;});

    // Look at the path in the same direction
    std::for_each(forward_path.cbegin(), forward_path.cend(),
            [&l,&count_forward](const sgNodeID_t &n){l.find(n) != l.cend() ? ++count_forward:count_forward;});


    // Select the path that has the most amount of matches
    if (count_reversed > count_forward) {
        oriented_read_paths[rid] = reversed_path;
    } else {
        oriented_read_paths[rid] = forward_path;
    }
    std::cout << "Oriented read path: " << std::endl;
    std::copy(oriented_read_paths[rid].cbegin(), oriented_read_paths[rid].cend(), std::ostream_iterator<sgNodeID_t>(std::cout, ", "));
    std::cout << std::endl;

}

void HaplotypeConsensus::build_line_path() {
    std::vector<sgNodeID_t> line_path;
    line_path.emplace_back(backbone[0]);

    for (int gap_number=1; gap_number < backbone.size(); gap_number++) {
        auto n1=backbone[gap_number-1];
        auto n2=backbone[gap_number];
        std::cout << "\n\nPrinting paths between "<<n1 << ", " << " and "<< n2 << ":\n";
        std::map<std::vector<sgNodeID_t>, uint32_t> gap_paths;
        for (const auto &rid :long_reads_in_backbone) {
            const auto &p = oriented_read_paths[rid];
            auto pn1 = std::find(p.cbegin(), p.cend(), n1);
            auto pn2 = std::find(p.cbegin(), p.cend(), n2);
//            if (std::distance(p.cbegin(), pn2) < std::distance(p.cbegin(), pn1)) std::swap(pn1,pn2);
            if (pn1 != p.cend() and pn2 != p.cend()){
                ++gap_paths[std::vector<sgNodeID_t>(pn1, pn2)];
            }
        }
        std::multimap<uint32_t, std::vector<sgNodeID_t>> most_common = flip_map(gap_paths);

        auto filled(false);
        for (const auto &p: most_common) {
            std::cout << p.first << ": ";
            std::copy(p.second.cbegin(), p.second.cend(), std::ostream_iterator<sgNodeID_t>(std::cout, ", "));

            if (!filled and p.first>1 and std::find(p.second.cbegin(), p.second.cend(), 0) == p.second.cend()){
                filled = true;
                auto it = ++p.second.cbegin();
                line_path.insert(line_path.end(), it, p.second.cend());
                std::cout << " <-- WINNER" << std::endl;
            } else {
                std::cout << std::endl;
            }
        }
        if (!filled){
            line_path.emplace_back(0);
        }
        line_path.emplace_back(n2);
    }
    std::cout << "\n\nFull line_path: ";
    std::copy(line_path.cbegin(), line_path.cend(), std::ostream_iterator<sgNodeID_t>(std::cout, ", "));
    std::cout << std::endl;
    backbone_filled_path = line_path;
}
