//
// Created by Bernardo Clavijo (EI) on 2019-02-27.
//

#include "HaplotypeConsensus.hpp"
#include <sglib/utilities/omp_safe.hpp>
#include <sglib/utilities/most_common_helper.hpp>

std::string HaplotypeConsensus::consensus_sequence(int disconnected_distance, int min_distance) {
    std::string consensus;
    std::vector<sgNodeID_t > current_path;

    for (uint32_t ni = 0; ni < backbone_filled_path.size(); ++ni) {
        const auto &n = backbone_filled_path[ni];
        if (n!=0) {
            current_path.emplace_back(n);
        }
        if (n == 0 or ni==backbone_filled_path.size()-1) {
            consensus+=SequenceGraphPath(ldg.sg, current_path).get_sequence();
            current_path.clear();
        }
        if (n == 0 and ni < backbone_filled_path.size()-1){
            int dist=disconnected_distance;
            for (const auto &fwl: ldg.get_fw_links(backbone_filled_path[ni-1])) {
                if (fwl.dest == backbone_filled_path[ni+1]) dist=fwl.dist;
            }
            consensus += std::string(std::max(dist,min_distance), 'N');
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

void HaplotypeConsensus::orient_read_path(uint64_t rid) {
    std::set<sgNodeID_t > l(backbone.cbegin(), backbone.cend());

    const auto &original_path = ws.long_read_mappers[0].read_paths[rid];
//    std::cout << "Original read path " << rid << ": " << std::endl;
//    std::copy(original_path.cbegin(), original_path.cend(), std::ostream_iterator<sgNodeID_t>(std::cout, ", "));
//    std::cout << std::endl;

    std::vector<sgNodeID_t> reversed_path;
    uint32_t count_reversed(0);
    uint32_t count_forward(0);
    // Look at the path reversed and with opposite signs
    std::transform(original_path.crbegin(), original_path.crend(), std::back_inserter(reversed_path), [](const sgNodeID_t &n) {return -n;});
    std::for_each(reversed_path.cbegin(), reversed_path.cend(),
            [&l,&count_reversed](const sgNodeID_t &n){l.find(n) != l.cend() ? ++count_reversed:count_reversed;});

    // Look at the path in the same direction
    std::for_each(original_path.cbegin(), original_path.cend(),
            [&l,&count_forward](const sgNodeID_t &n){l.find(n) != l.cend() ? ++count_forward:count_forward;});


    // Select the path that has the most amount of matches
    if (count_reversed > count_forward) {
        oriented_read_paths[rid] = reversed_path;
    } else {
        oriented_read_paths[rid] = ws.long_read_mappers[0].read_paths[rid];
    }
//    std::cout << "Oriented read path " << rid << ": " << std::endl;
//    std::copy(oriented_read_paths[rid].cbegin(), oriented_read_paths[rid].cend(), std::ostream_iterator<sgNodeID_t>(std::cout, ", "));
//    std::cout << std::endl;

}

void HaplotypeConsensus::build_line_path(int min_votes, int min_path_nodes) {
    std::vector<sgNodeID_t> final_line_path;
    final_line_path.emplace_back(backbone[0]);
    std::vector<std::vector<sgNodeID_t>> thread_line_paths(backbone.size()-1);
#pragma omp parallel
    { // Parallel section for all gaps
#pragma omp for
        for (int gap_number = 1; gap_number < backbone.size(); gap_number++) {
            auto n1 = backbone[gap_number - 1];
            auto n2 = backbone[gap_number];
            std::cout << "\n\nPrinting paths between " << n1 << ", " << " and " << n2 << ":\n";
            std::map<std::vector<sgNodeID_t>, uint32_t> gap_paths;
            auto read_path = oriented_read_paths.cbegin();
            for (int i = 0; read_path != oriented_read_paths.end(); i++, ++read_path) {
                const auto &p = *read_path;
                uint32_t pn1 = p.size();
                uint32_t pn2 = p.size();
                for (uint32_t path_node_pos = 0; path_node_pos < p.size(); path_node_pos++) {
                    if (p[path_node_pos] == n1 and pn1 == p.size()) {
                        pn1 = path_node_pos;
                    }
                    if (p[path_node_pos] == n2 and pn2 == p.size()) {
                        pn2 = path_node_pos;
                    }
                }
                if (pn1 > pn2) {
//                std::cout << n1 << " and " << n2 << " nodes from backbone are at " << std::distance(p.cbegin(), pn1) << " and " << std::distance(p.cbegin(), pn2) << " from start, they are wrongly oriented in paths!" << std::endl;
                    continue; // TODO: Some paths seem to be wrongly oriented?!
                }
                if (pn1 != p.size() and pn2 != p.size()) {
                    ++gap_paths[std::vector<sgNodeID_t>(p.begin()+pn1, p.begin()+pn2)];
                }
            }
            std::multimap<uint32_t, std::vector<sgNodeID_t>> most_common = flip_map(gap_paths);

            auto filled(false);
            for (auto p = most_common.crbegin(); p != most_common.crend(); ++p) {
                std::cout << p->first << ": ";
                std::copy(p->second.cbegin(), p->second.cend(), std::ostream_iterator<sgNodeID_t>(std::cout, ", "));

                // If there's a path with more than min_votes vote and it doesn't contain a gap, then this path is the winner
                if (!filled and p->first > min_votes and
                    std::find(p->second.cbegin(), p->second.cend(), 0) == p->second.cend()) {
                    filled = true;
                    auto it = ++p->second.cbegin();
                    thread_line_paths[gap_number].insert(thread_line_paths[gap_number].end(), it, p->second.cend());
                    std::cout << " <-- WINNER" << std::endl;
                } else {
                    std::cout << std::endl;
                }
            }

            // If no path without gaps and more than min_votes was found then collect all the "common" nodes as partial solutions
            // and place a gap in between the uncommon nodes
            if (!filled) {
                std::vector<std::vector<sgNodeID_t>> winners;
                for (const auto &mcp : most_common) {
                    // Only consider paths with more than min_path_nodes nodes for partial solutions, typically
                    // the first shared node is the first node in the backbone and the second is a gap, this if
                    // avoids those cases but will need to be changed if the path doesn't include the first backbone node
                    if (mcp.second.size() > min_path_nodes) {
                        winners.emplace_back(mcp.second);
                    }
                }
                if (winners.empty()) {
                    thread_line_paths[gap_number].emplace_back(0);
                    thread_line_paths[gap_number].emplace_back(n2);
                    continue;
                }
                std::sort(winners.begin(), winners.end(),
                          [](const std::vector<sgNodeID_t> &a, const std::vector<sgNodeID_t> &b) {
                              return a.size() < b.size();
                          });
                std::vector<sgNodeID_t> shared_winner_nodes;
                for (int pos = 1; pos < winners[0].size() and winners[0][pos] != 0; ++pos) {
                    bool sharing = true;
                    auto &current_node = winners[0][pos];
                    for (const auto &w : winners) {
                        if (w[pos] != current_node or w[pos] == 0) {
                            sharing = false;
                            break;
                        }
                    }
                    if (sharing) shared_winner_nodes.emplace_back(current_node);
                    else break;
                }

                std::vector<sgNodeID_t> back_shared_winner_nodes;
                for (int pos = 0; pos < winners[0].size() and winners[0][winners[0].size() - pos - 1] != 0; ++pos) {
                    bool sharing = true;
                    auto &current_node = winners[0][winners[0].size() - pos - 1];
                    for (const auto &w : winners) {
                        if (w[w.size() - pos - 1] != current_node or w[w.size() - pos - 1] == 0) {
                            sharing = false;
                            break;
                        }
                    }
                    if (sharing) back_shared_winner_nodes.emplace_back(current_node);
                    else break;
                }

                thread_line_paths[gap_number].insert(thread_line_paths[gap_number].end(), shared_winner_nodes.begin(), shared_winner_nodes.end());
                thread_line_paths[gap_number].emplace_back(0);  // TODO: Check why there are more N's now with same number of gaps!
                thread_line_paths[gap_number].insert(thread_line_paths[gap_number].end(), back_shared_winner_nodes.rbegin(), back_shared_winner_nodes.rend());

//            std::cout << "Partial line path from shared nodes in paths: " << std::endl;
//            std::cout << "First shared winners: "; std::copy(shared_winner_nodes.cbegin(), shared_winner_nodes.cend(), std::ostream_iterator<sgNodeID_t>(std::cout, ", "));
//            std::cout << std::endl;
//            std::cout << "0, ";
//            std::cout << std::endl;
//            std::cout << "Last shared winners: "; std::copy(back_shared_winner_nodes.crbegin(), back_shared_winner_nodes.crend(), std::ostream_iterator<sgNodeID_t>(std::cout, ", "));
//            std::cout << std::endl;

            }
            thread_line_paths[gap_number].emplace_back(n2);
        }

#pragma omp single
        {
            for (uint32_t tid=0; tid < backbone.size()-1; tid++) {
                final_line_path.insert(final_line_path.end(), thread_line_paths[tid].begin(), thread_line_paths[tid].end());
            }
        }
    }
    std::cout << "\n\nFull line_path: ";
    std::copy(final_line_path.cbegin(), final_line_path.cend(), std::ostream_iterator<sgNodeID_t>(std::cout, ", "));
    std::cout << std::endl;
    backbone_filled_path = final_line_path;
}
