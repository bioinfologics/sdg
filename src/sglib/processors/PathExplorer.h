//
// Created by Ben Ward (EI) on 22/03/2018.
//

#ifndef BSG_PATHEXPLORER_H
#define BSG_PATHEXPLORER_H

#include <sglib/graph/SequenceGraph.hpp>

class PathExplorer {
public:
    explicit PathExplorer(SequenceGraph& sg) : sg(sg) {}
    std::vector<SequenceGraphPath> collect_paths(const sgNodeID_t& seed, const sgNodeID_t& target,
                                                 unsigned int size_limit = 0, unsigned int edge_limit = 0) const;

    std::vector<SequenceGraphPath> collect_paths(const sgNodeID_t &seed, const sgNodeID_t &target,
                                                 const std::string& query, unsigned int flank = 500) const;

    int find_best_path(SequenceGraphPath& result,
                       const sgNodeID_t& seed,
                       const sgNodeID_t& target,
                       const std::string& query,
                       unsigned int flank = 500) const;

private:
    SequenceGraph& sg;
};


#endif //BSG_PATHEXPLORER_H
