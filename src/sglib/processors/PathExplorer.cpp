//
// Created by Ben Ward (EI) on 22/03/2018.
//

#include <stack>
#include <sglib/graph/SequenceGraphPath.hpp>
#include <sglib/processors/PathExplorer.h>
#include <sglib/aligners/algorithms/SmithWaterman.h>

std::vector<SequenceGraphPath> PathExplorer::collect_paths(const sgNodeID_t &seed, const sgNodeID_t &target,
                                                           const std::string& query, const unsigned int flank) const {
    const unsigned int size_limit = (unsigned int) query.size() + flank;
    return collect_paths(seed, target, size_limit);
}

std::vector<SequenceGraphPath> PathExplorer::collect_paths(const sgNodeID_t& seed,
                                                           const sgNodeID_t& target,
                                                           const unsigned int size_limit,
                                                           const unsigned int edge_limit) const {

    struct visitor {
        sgNodeID_t node;
        uint dist;
        uint path_length;
        visitor(sgNodeID_t n, uint d, uint p) : node(n), dist(d), path_length(p) {}
    };

    std::vector<SequenceGraphPath> collected_paths;
    std::vector<sgNodeID_t> current_path;

    std::stack<visitor> to_visit;
    to_visit.emplace(seed, 0, 0);
    std::set<sgNodeID_t> visited;

    while (!to_visit.empty()) {
        const auto activeNode(to_visit.top());
        to_visit.pop();
        if (visited.find(activeNode.node) == visited.end() and (activeNode.path_length < edge_limit or edge_limit==0) and (activeNode.dist < size_limit or size_limit==0) )
        {
            visited.emplace(activeNode.node);

            if (current_path.size() > activeNode.path_length) {
                current_path[activeNode.path_length] = activeNode.node;
                current_path.resize(activeNode.path_length + 1);
            } else {
                current_path.emplace_back(activeNode.node);
            }

            for (const auto &l : sg.get_fw_links(activeNode.node)) {
                // For each child of the current node, create a child visitor object, and enter it into the parent map.
                visitor child { l.dest, activeNode.dist + uint(sg.nodes[std::abs(l.dest)].sequence.length()), activeNode.path_length + 1 };

                // If I've found the target along this dfs path, I need to collect it and put it into a SequenceGraphPath object.
                if (child.node == target) {
                    current_path.emplace_back(target);
                    collected_paths.emplace_back(sg, current_path);
                } else {
                    to_visit.emplace(child);
                }
            }
        }
    }
    return collected_paths;
}

int PathExplorer::find_best_path(SequenceGraphPath& result,
                                 const sgNodeID_t& seed,
                                 const sgNodeID_t& target,
                                 const std::string& query,
                                 const unsigned int flank) const {

    const auto paths { collect_paths(seed, target, query, flank) };

    if (!paths.empty()) {
        if (paths.size() == 1) {
            result = paths[0];
            return 0;
        } else {
            // Determine the best path of several, using the Smith Waterman alignment.
            sglib::alignment::algorithms::SmithWaterman<int> sw { paths.front().get_sequence().size(), query.size() };
            int maxscore { 0 };
            std::vector<SequenceGraphPath>::const_iterator bestpath;
            for (auto path {paths.cbegin()}; path != paths.cend(); ++path) {
                const auto alignres { sw.run(path -> get_sequence(), query, -10, -1) };
                if (std::get<0>(alignres) > maxscore) {
                    maxscore = std::get<0>(alignres);
                    bestpath = path;
                }
            }
            result = *bestpath;
        }
    } else {
        return 1;
    }

}