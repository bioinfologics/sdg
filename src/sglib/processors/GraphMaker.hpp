//
// Created by Bernardo Clavijo (EI) on 05/07/2018.
//

// This class IS a re-invention of the wheel. It creates a DBG from a number of different sources and puts it into a
// SequenceGraph.


#ifndef BSG_GRAPHMAKER_HPP
#define BSG_GRAPHMAKER_HPP


#include <sglib/SequenceGraph.hpp>

/**
 * This kmer class can give possible neighbours FW and BW to help build a graph more easily.
 *
 */
class GraphMakerKmer{
public:
    GraphMakerKmer(uint64_t _kmervalue, uint8_t _k):kmervalue(_kmervalue),k(_k){};

    uint64_t kmervalue;
    uint8_t k;
};

class GraphMaker {
public:
    GraphMaker(SequenceGraph & _sg): sg(_sg){};
    void new_graph_from_kmerset_trivial(const std::unordered_set<uint64_t> & kmerset,uint8_t k);
private:
    SequenceGraph & sg;
};


#endif //BSG_GRAPHMAKER_HPP
