//
// Created by Luis Yanes (EI) on 22/03/2018.
//

#ifndef BSG_SEQUENCEGRAPHPATH_HPP
#define BSG_SEQUENCEGRAPHPATH_HPP


#include <sglib/graph/SequenceGraph.hpp>

class SequenceGraph;
class SequenceGraphPath {
public:
    std::vector<sgNodeID_t> nodes;
    explicit SequenceGraphPath(SequenceGraph & _sg, const std::vector<sgNodeID_t> _nodes={})  : sg(_sg) ,nodes(_nodes) {};

    SequenceGraphPath(const SequenceGraphPath& sgp) : nodes(sgp.nodes), sg(sgp.sg) {};

    SequenceGraphPath& operator=(const SequenceGraphPath &other);

    std::string get_fasta_header(bool use_oldnames = false) const;
    std::string get_sequence() const;
    std::vector<Link> get_next_links();
    void reverse();
    bool is_canonical();
    std::set<sgNodeID_t> make_set_of_nodes() const;
    bool operator==(const SequenceGraphPath& rhs) const;
    bool operator<(const SequenceGraphPath& rhs) const;
    bool append_to_path(sgNodeID_t newnode);
    bool extend_if_coherent(SequenceGraphPath s);
    void clear() {
        nodes.clear();
    };

private:
    const SequenceGraph& sg;
};


#endif //BSG_SEQUENCEGRAPHPATH_HPP
