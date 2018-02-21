//
// Created by Bernardo Clavijo (EI) on 10/11/2017.
//

#ifndef SG_SCAFFOLDER_HPP
#define SG_SCAFFOLDER_HPP

#include <sglib/mappers/LinkedReadMapper.hpp>
#include "SequenceGraph.h"
#include "PairedReadMapper.h"
#include "KmerCompressionIndex.hpp"

class Scaffolder {

public:
    /**
     * @brief
     * The Scaffolder object provides functions to locate or resolve regions
     * of repetitions and haplotype originated bubbles.
     * It requires a valid SequenceGraph, a vector of ReadMappers and a KmerCompressionIndex.
     * @param _sg A SequenceGraph to scaffold
     * @param _rms A vector of PairedReadmapper, containing the mapping of reads to _sg
     * @param _kci A KMerCompressionIndex object
     */
    Scaffolder(SequenceGraph &_sg, std::vector<PairedReadMapper> & _rms,  std::vector<LinkedReadMapper> & _lrms, KmerCompressionIndex &_kci) : sg(_sg),rmappers(_rms),lrmappers(_lrms),kci(_kci){};

    void pop_unsupported_shortbubbles();
    void expand_bubbly_subgraphs();
    std::vector<SequenceSubGraph> get_all_bubbly_subgraphs(uint32_t maxsubgraphs=0);
    std::vector<std::pair<sgNodeID_t,sgNodeID_t>> get_all_haplotype_pairs(uint32_t maxpairs=0);

    void find_canonical_repeats();

    ///
    /// \param source
    /// \param dest
    /// \return count of all the reads from all mappers, linking source -> dest or source -> -dest
    uint64_t count_reads_linking(sgNodeID_t source, sgNodeID_t dest);

    ///
    /// \param source
    /// \param lib the library (i.e. mapper index)
    /// \return a vecttor with dests (signed) and how many reads link to them
    std::vector<std::pair<sgNodeID_t,uint64_t>> all_read_links(sgNodeID_t source, unsigned lib);

    //2: find unsatisfied connections
    //find read support breakpoints
    //3: path finder?


    //This is a trivial strategy for the incorporation of long reads:
    //1) create a path for each long read
    //2) combine paths (OL)

    void findLongreadPaths();
    void joinLongreadPaths();
    void applyLonreadPaths();

    SequenceGraph &sg;
    std::vector<PairedReadMapper> &rmappers;
    std::vector<LinkedReadMapper> &lrmappers;
    KmerCompressionIndex &kci;



};


#endif //SG_SCAFFOLDER_HPP
