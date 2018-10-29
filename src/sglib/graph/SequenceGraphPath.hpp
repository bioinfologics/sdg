//
// Created by Luis Yanes (EI) on 22/03/2018.
//

#ifndef BSG_SEQUENCEGRAPHPATH_HPP
#define BSG_SEQUENCEGRAPHPATH_HPP

#include <string>
#include <vector>
#include <set>
#include <sglib/types/GenericTypes.hpp>

class SequenceSubGraph;
class SequenceGraph;
class SequenceGraphPath {
public:
    explicit SequenceGraphPath(SequenceGraph & _sg, const std::vector<sgNodeID_t> _nodes={})  : sg(_sg) ,nodes(_nodes) {};
    explicit SequenceGraphPath(const SequenceGraph &_sg, const std::vector<sgNodeID_t> _nodes={}) : sg(_sg), nodes(_nodes) {};

    SequenceGraphPath(const SequenceGraphPath& sgp) : nodes(sgp.nodes), sg(sgp.sg) {};

    SequenceGraphPath& operator=(const SequenceGraphPath &other);

    std::string get_fasta_header(bool use_oldnames = false) const;
    std::string get_sequence() const;
    size_t get_sequence_size_fast();
    std::vector<Link> get_next_links();
    void reverse();
    bool is_canonical();
    std::set<sgNodeID_t> make_set_of_nodes() const;
    bool operator==(const SequenceGraphPath& rhs) const;
    bool operator<(const SequenceGraphPath& rhs) const;
    bool append_to_path(sgNodeID_t newnode);
    bool extend_if_coherent(SequenceGraphPath s){};
    void clear() {
        nodes.clear();
    };

    bool is_unitig();
    std::vector<sgNodeID_t >& getNodes() {return nodes;}
    const std::vector<sgNodeID_t >& getNodes() const {return nodes;}

    std::vector<sgNodeID_t> nodes;

private:
    const SequenceGraph& sg;
};


#endif //BSG_SEQUENCEGRAPHPATH_HPP
