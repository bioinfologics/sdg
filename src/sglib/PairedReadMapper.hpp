//
// Created by Bernardo Clavijo (EI) on 03/11/2017.
//

#ifndef SG_PAIREDREADMAPPER_HPP
#define SG_PAIREDREADMAPPER_HPP


#include "SequenceGraph.hpp"
enum prmReadType {prmPE, prmLMP, prm10x};
const std::string prmReadTypeDesc[]={"Paired End", "Long Mate Pair", "10x Linked Reads"};

class ReadMapping {
public:
    uint32_t read_id;
    sgNodeID_t node;
    int32_t first_pos;
    int32_t last_pos;
    int32_t unique_matches;
};

class PairedReadMapper {
public:
    PairedReadMapper(SequenceGraph &_sg) : sg(_sg){
        reads_in_node.resize(sg.nodes.size());
    };
    void map_reads(std::string filename1, std::string filename2, prmReadType read_type=prmPE);

    SequenceGraph & sg;
    std::vector<std::vector<ReadMapping>> reads_in_node;
};


#endif //SG_PAIREDREADMAPPER_HPP
