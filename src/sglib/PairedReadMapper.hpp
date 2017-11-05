//
// Created by Bernardo Clavijo (EI) on 03/11/2017.
//

#ifndef SG_PAIREDREADMAPPER_HPP
#define SG_PAIREDREADMAPPER_HPP


#include "SequenceGraph.hpp"
#include "KMerIDX.h"
#include "FileReader.h"
#include "SMR.h"

enum prmReadType {prmPE, prmLMP, prm10x};
const std::string prmReadTypeDesc[]={"Paired End", "Long Mate Pair", "10x Linked Reads"};


class ReadMapping {
public:
    sgNodeID_t node;
    uint32_t read_id;
    int32_t first_pos;
    int32_t last_pos;
    int32_t unique_matches;
    bool operator==(const ReadMapping &other){
        return this==&other;
    };
    bool operator<(const ReadMapping &other) const {
        return node<other.node;
    };
    void merge(const ReadMapping &other){};
};

//class NodeReadContainer {
//    sgNodeID_t node;
//    std::vector<ReadMapping> reads;
//
//    NodeReadContainer(sgNodeID_t n=0): node(n) {}
//
//    const bool operator<(const NodeReadContainer& other) const {
//        return (node < other.node);
//    }
//
//    const bool operator>(const NodeReadContainer& other) const {
//        return (node > other.node);
//    }
//
//    const bool operator==(const NodeReadContainer &other) const {
//        return ( node == other.node );
//    }
//
//    void merge(const NodeReadContainer &other) {
//        reads.reserve(reads.size()+other.reads.size());
//        reads.insert(reads.end(), other.reads.begin(), other.reads.end());
//    }
//    NodeReadContainer max() { return {}; }
//
//    friend std::ostream& operator<<(std::ostream& os, const NodeReadContainer& mr){
//        os << link.tag << "\t" << link.contig << ":" << (int) link.count;
//        return os;
//    }
//    friend std::istream& operator>>(std::istream& is, const NodeReadContainer& mr) {
//        is.read((char*)&link, sizeof(link));
//        return is;
//    }
//};

class PairedReadMapper {
public:
    PairedReadMapper(SequenceGraph &_sg) : sg(_sg){
        reads_in_node.resize(sg.nodes.size());
    };
    void map_reads(std::string filename1, std::string filename2, prmReadType read_type=prmPE);
    uint64_t process_reads_from_file(uint8_t k, uint16_t min_matches, std::vector<KmerIDX> &unique_kmers, std::string filename, uint64_t offset);

    SequenceGraph & sg;
    std::vector<std::vector<ReadMapping>> reads_in_node;
    std::vector<sgNodeID_t> read_to_node;
};


#endif //SG_PAIREDREADMAPPER_HPP
