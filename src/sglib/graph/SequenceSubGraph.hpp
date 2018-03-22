//
// Created by Luis Yanes (EI) on 22/03/2018.
//

#ifndef BSG_SEQUENCESUBGRAPH_HPP
#define BSG_SEQUENCESUBGRAPH_HPP

#include <sglib/graph/SequenceGraph.h>

class SequenceGraph;
class SequenceGraphPath;
class SequenceSubGraph {
public:
    std::vector<sgNodeID_t> nodes;
    explicit SequenceSubGraph(SequenceGraph & _sg, std::vector<sgNodeID_t> _nodes={})  : sg(_sg) ,nodes(_nodes) {};
    SequenceGraphPath make_path(); // Returns empty path if not linear.

    void write_to_gfa(std::string filename);
private:
    SequenceGraph& sg;
};


#endif //BSG_SEQUENCESUBGRAPH_HPP
