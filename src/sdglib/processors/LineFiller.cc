//
// Created by Gonzalo Garcia (EI) on 2020-07-07.
//

#include <sdglib/workspace/WorkSpace.hpp>
#include "LineFiller.h"

uint32_t LineFiller::score_function(std::vector<sgNodeID_t> path) {

    // Get set of absolute nodes
    std::vector<sgNodeID_t> abs_nodes;
    abs_nodes.reserve(path.size());
    for(const auto &node: path){
        abs_nodes.push_back(abs(node));
    }

    uint32_t path_score = 0;
    for (const auto &node: path){
        path_score+=node_matches[abs(node)];
    }
    return path_score;
}

void LineFiller::populate_matches(){
    sdglib::OutputLog() << "Starting "<< std::endl;
    node_matches.resize(ws.sdg.nodes.size(), 0);
    sdglib::OutputLog() << "Vector resized "<< std::endl;
    for (const auto &matches: lrr.read_perfect_matches){
        for (const auto& match: matches){
            node_matches[abs(match.node)]+=match.size;
        }
    }
    sdglib::OutputLog() << "Done "<< std::endl;
}

std::vector<sgNodeID_t> LineFiller::line_fill(std::vector<sgNodeID_t> anchor_path){
    std::vector<sgNodeID_t> final_path;
    for (auto i=0; i<anchor_path.size()-1; ++i){
        auto n1 = anchor_path[i];
        auto n2 = anchor_path[i+1];
        std::cout << "Filling gap between "<< n1 << " and " << n2 << std::endl;
        final_path.push_back(n1);

        auto paths_between = ws.sdg.find_all_paths_between(n1, n2, 10000, 100, false);
        std::cout << "Paths between "<< n1 << " and " << n2 << " --> "<< paths_between.size() << std::endl;

        if (paths_between.size() == 0){
            continue;
        }
        else if (paths_between.size() == 1){
            final_path.insert(final_path.end(), paths_between[0].nodes.begin(), paths_between[0].nodes.end());
            continue;
        }
        else{
            auto max_score = 0;
            auto max_score_index = 0;
            for (auto p=0; p<paths_between.size(); ++i){
                auto score = score_function(paths_between[p].nodes);
                if (score>max_score){
                    max_score = score;
                    max_score_index = p;
                }
            }
            final_path.insert(final_path.end(), paths_between[max_score_index].nodes.begin(), paths_between[max_score_index].nodes.end());
        }

    }
    final_path.push_back(anchor_path[anchor_path.size()]);
    return final_path;
}