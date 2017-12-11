//
// Created by Bernardo Clavijo (EI) on 03/11/2017.
//

#ifndef SG_PAIREDREADMAPPER_HPP
#define SG_PAIREDREADMAPPER_HPP

#include <map>

#include "SequenceGraph.hpp"
#include "sglib/factories/KMerIDXFactory.h"
#include "sglib/readers/FileReader.h"
#include "sglib/readers/SequenceGraphReader.h"
#include "SMR.h"

enum prmReadType {prmPE, prmLMP, prm10x};
const std::string prmReadTypeDesc[]={"Paired End", "Long Mate Pair", "10x Linked Reads"};
typedef uint32_t prm10xTag_t;

class ReadMapping {
public:
    ReadMapping(){
        //just clean the structure, OSX doesn't give you clean memory
        bzero(this, sizeof(ReadMapping));
    }
    bool operator==(const ReadMapping &other){
        return this==&other;
    };
    bool operator<(const ReadMapping &other) const {
        if (node!=other.node) return node<other.node;
        return read_id<other.read_id;
    };
    void merge(const ReadMapping &other){};

    sgNodeID_t node;
    uint64_t read_id;
    int32_t first_pos;
    int32_t last_pos;
    int32_t unique_matches;
    bool rev=false;

};

/**
 * @brief A mapper for different kinds of paired short reads (Illumina).
 *
 * The mapper supports so far PE and 10x reads.
 * 10x reads must be trimmed and have the tags at the end of the read name.
 */
class PairedReadMapper {
public:
    PairedReadMapper(SequenceGraph &_sg) : sg(_sg){
        std::cout << "_sg size " << _sg.nodes.size();
        std::cout << "sg size " << sg.nodes.size();
        reads_in_node.resize(sg.nodes.size());
        std::cout << "reads_in_node size; " << reads_in_node.size() << std::endl;
    };
    void map_reads(std::string , std::string , prmReadType , uint64_t );
    void remove_obsolete_mappings();
    void remap_reads();
    uint64_t process_reads_from_file(uint8_t, uint16_t, std::vector<KmerIDX> &, std::string , uint64_t, bool );
    void save_to_disk(std::string filename);
    void load_from_disk(std::string filename);
    void print_stats();

    SequenceGraph & sg;
    std::string read1filename,read2filename;
    prmReadType readType;
    uint64_t memlimit;
    std::vector<std::vector<ReadMapping>> reads_in_node;
    std::vector<sgNodeID_t> read_to_node;//id of the main node if mapped, set to 0 to remap on next process
    std::vector<prm10xTag_t> read_to_tag;
};


#endif //SG_PAIREDREADMAPPER_HPP
