//
// Created by Luis Yanes (EI) on 22/03/2018.
//

#ifndef BSG_SEQUENCESUBGRAPH_HPP
#define BSG_SEQUENCESUBGRAPH_HPP

#include "sglib/graph/SequenceGraph.hpp"
#include "sglib/graph/SequenceGraphPath.hpp"
class SequenceGraph;

class SequenceSubGraph {
public:
    std::vector<sgNodeID_t> nodes;
    explicit SequenceSubGraph(const SequenceGraph & _sg, std::vector<sgNodeID_t> _nodes={})  : sg(_sg) ,nodes(_nodes) {};
    explicit SequenceSubGraph(const SequenceGraph & _sg, std::vector<nodeVisitor> nodes) : sg(_sg) {
        for (const auto &nv : nodes) {
            this->nodes.emplace_back(nv.node);
        }
    };
//    SequenceGraphPath make_path(){}; // Returns empty path if not linear.

    void write_to_gfa(std::string filename);
    uint64_t total_size() const ;

private:
    const SequenceGraph& sg;
};


#endif //BSG_SEQUENCESUBGRAPH_HPP
