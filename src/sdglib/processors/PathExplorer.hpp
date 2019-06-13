//
// Created by Ben Ward (EI) on 22/03/2018.
//

#ifndef BSG_PATHEXPLORER_H
#define BSG_PATHEXPLORER_H

#include <sdglib/graph/SequenceDistanceGraph.hpp>

class PathExplorer {
public:
    explicit PathExplorer(const SequenceDistanceGraph &sg) : sg(sg) {}
    std::vector<SequenceGraphPath> collect_paths(sgNodeID_t seed, sgNodeID_t target, unsigned int size_limit = 0, unsigned int edge_limit = 0) const;
    std::vector<SequenceGraphPath> collect_paths(sgNodeID_t seed, sgNodeID_t target, const std::string& query, unsigned int flank = 500) const;
    bool find_best_path(SequenceGraphPath& result, sgNodeID_t seed, sgNodeID_t target, const std::string& query, unsigned int flank = 500) const;

private:
    const SequenceDistanceGraph& sg;
};


#endif //BSG_PATHEXPLORER_H
