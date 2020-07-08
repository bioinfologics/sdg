//
// Created by Gonzalo Garcia (EI) on 2020-07-07.
//

#include <sdglib/workspace/WorkSpace.hpp>
#include "LineFiller.h"

std::vector<sgNodeID_t> LineFiller::score_function(sgNodeID_t n1, sgNodeID_t n2, std::vector<SequenceDistanceGraphPath> paths) {

    // Find all reads spanning both nodes

    auto reads_n1 = lrr.node_reads[abs(n1)];
    auto reads_n2 = lrr.node_reads[abs(n2)];
//    std::cout << "Reads in n1: "<< reads_n1.size()<< " reads in n2: "<< reads_n2.size() << std::endl;

    // Get all matches for both ends of the gap
    std::vector<int64_t> shared_reads;
    std::set_intersection(reads_n1.begin(), reads_n1.end(), reads_n2.begin(), reads_n2.end(),  std::back_inserter(shared_reads), [](int64_t r1, int64_t r2){return abs(r1) < abs(r2);});
//    std::cout << "Shared reads: " << shared_reads.size() << std::endl;

    std::vector<PerfectMatch> shared_matches;
    for (const auto& rid: shared_reads){
        for (const auto& match: lrr.read_perfect_matches[abs(rid)]){
            shared_matches.push_back(match);
        }
    }
//    std::cout << "Shared matches: "<< shared_matches.size() << std::endl;

    int32_t max_score = 0;
    int32_t max_score_index = 0;
    int pos = 0;
    for (const auto& path: paths){
        // Absolute node collection for comparison
        int32_t score = 0;
        std::vector<sgNodeID_t> abs_nodes;
        for (const auto& node: path.nodes){
            abs_nodes.push_back(abs(node));
        }

        for (const auto& match: shared_matches){
            if (std::find(abs_nodes.begin(), abs_nodes.end(), abs(match.node)) != abs_nodes.end()){
                score+=match.size;
            }
        }

        if (score>max_score){
            max_score = score;
            max_score_index = pos;
        }

        pos++;
    }
//    std::cout << "Max score: "<< max_score<< " max position: "<< max_score_index << std::endl;
    return paths[max_score_index].nodes;
}

std::vector<sgNodeID_t> LineFiller::line_fill(std::vector<sgNodeID_t> anchor_path){
    std::vector<sgNodeID_t> final_path;
    for (auto i=0; i<anchor_path.size()-1; ++i){
        auto n1 = anchor_path[i];
        auto n2 = anchor_path[i+1];
//        std::cout << "Filling gap between "<< n1 << " and " << n2 << std::endl;
        final_path.push_back(n1);

        auto paths_between = ws.sdg.find_all_paths_between(n1, n2, 10000, 100, false);
//        std::cout << "Paths between "<< n1 << " and " << n2 << " --> "<< paths_between.size() << std::endl;

        if (paths_between.size() == 0){
            continue;
        }
        else if (paths_between.size() == 1){
            final_path.insert(final_path.end(), paths_between[0].nodes.begin(), paths_between[0].nodes.end());
            continue;
        }
        else{
            // More than 1 path insert the best scoring
            auto winner_path = score_function(n1, n2, paths_between);
            final_path.insert(final_path.end(), winner_path.begin(), winner_path.end());
        }

    }
    //final_path.push_back(anchor_path[anchor_path.size()]);
    return final_path;
}

std::vector<std::vector<sgNodeID_t >> LineFiller::fill_all_paths(std::vector<std::vector<sgNodeID_t >> lines){
    std::vector<std::vector<sgNodeID_t>> final_lines(lines.size());
    auto total_ready = 0;
#pragma omp parallel for
    for (auto l=0; l<lines.size(); ++l){
        auto line = lines[l];
        std::cout << "Processing line " << total_ready << "/" << lines.size() <<std::endl;
        final_lines[l] = line_fill(line);
        total_ready++;
    }
    return final_lines;
};