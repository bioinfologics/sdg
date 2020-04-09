//
// Created by Bernardo Clavijo (EI) on 07/04/2020.
//

#pragma once


#include <sdglib/workspace/WorkSpace.hpp>
#include <sdglib/mappers/LongReadsRecruiter.hpp>

/*
 * This tries to choose the right path between two haplotype specific nodes, by using whatever data we can put to it.
 */
enum PFSEType {PFLongRead, PFshortsingle, PFshortpair};
class PFSequenceEvidence{
public:
    PFSequenceEvidence(std::string _seq, PFSEType _evtype, uint8_t _dsidx, uint8_t _rid):
        seq(_seq),evtype(_evtype),dsidx(_dsidx),rid(_rid){};
    std::string seq;
    PFSEType evtype;
    uint8_t dsidx;
    uint64_t rid;
};

class PathFinder {
public:
    PathFinder(WorkSpace &_ws, sgNodeID_t _n1, sgNodeID_t _n2, uint8_t _k):
        ws(_ws),n1(_n1),n2(_n2),k(_k){

    };

    //TODO: recive a MLDG and a LRR and chop/rc reads to extract sequence between n1 and n2
    void load_lrseqs(DistanceGraph &dg, const LongReadsRecruiter & lrr, int ovl_extension=300);
    std::string lrseqs_as_fasta();
    //Strategy #1: index the paths, then score each read against the paths to get the winner.
    void index_paths(std::vector<SequenceDistanceGraphPath> _paths);
    std::vector<std::vector<uint64_t>> seq_to_pathpos(uint16_t path_id, std::string seq);

    //Strategy #2: index the reads, advance through paths using a greedy strategy supported by reads
    void index_seqs();

    SequenceDistanceGraphPath find_path_with_lrseqs();

    //TODO: Strategy #3: local TAP
    std::unordered_map<uint64_t,std::vector<std::pair<uint16_t ,uint64_t >>> kmerpos; //kmer -> [(path_id, pos)]
    WorkSpace & ws;
    sgNodeID_t n1,n2;
    std::vector<SequenceDistanceGraphPath> paths;
    uint8_t k;
    std::vector<PFSequenceEvidence> seqs;
};


