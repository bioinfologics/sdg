//
// Created by Luis Yanes (EI) on 22/03/2018.
//

#ifndef BSG_SEQUENCEGRAPHPATH_HPP
#define BSG_SEQUENCEGRAPHPATH_HPP

#include <string>
#include <vector>
#include <set>
#include <sdglib/types/GenericTypes.hpp>

class SequenceSubGraph;
class SequenceDistanceGraph;
class SequenceDistanceGraphPath {
public:
    explicit SequenceDistanceGraphPath(SequenceDistanceGraph & _sg, const std::vector<sgNodeID_t> _nodes={})  : sg(_sg) ,nodes(_nodes) {};
    explicit SequenceDistanceGraphPath(const SequenceDistanceGraph &_sg, const std::vector<sgNodeID_t> _nodes={}) : sg(_sg), nodes(_nodes) {};

    SequenceDistanceGraphPath(const SequenceDistanceGraphPath& sgp) : nodes(sgp.nodes), sg(sgp.sg) {};

    SequenceDistanceGraphPath& operator=(const SequenceDistanceGraphPath &other);

    friend std::ostream& operator<<(std::ostream &os, const SequenceDistanceGraphPath &sdgp);

    std::string get_fasta_header(bool use_oldnames = false) const;
    std::string sequence() const;
    size_t get_sequence_size_fast() const;
    std::vector<Link> get_next_links();
    void reverse();
    bool is_canonical();
    std::set<sgNodeID_t> make_set_of_nodes() const;
    bool operator==(const SequenceDistanceGraphPath& rhs) const;
    bool operator<(const SequenceDistanceGraphPath& rhs) const;

    bool extend_if_coherent(SequenceDistanceGraphPath s){};
    void clear() {
        nodes.clear();
    };

    bool is_valid(bool allow_ns=true);
    bool is_unitig();
    std::vector<sgNodeID_t >& getNodes() {return nodes;}
    const std::vector<sgNodeID_t >& getNodes() const {return nodes;}

    std::vector<sgNodeID_t> nodes;

    const SequenceDistanceGraph& sg;
};


#endif //BSG_SEQUENCEGRAPHPATH_HPP
