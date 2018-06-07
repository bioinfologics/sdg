//
// Created by Bernardo Clavijo (EI) on 12/05/2018.
//

#ifndef BSG_PAIREDREADMAPPER_HPP
#define BSG_PAIREDREADMAPPER_HPP

#include <map>

#include "sglib/types/MappingTypes.hpp"
#include "sglib/factories/KMerIDXFactory.h"
#include "sglib/readers/SequenceGraphReader.h"
#include "sglib/SMR.h"
#include <sglib/datastores/PairedReadsDatastore.hpp>

class PairedReadConnectivityDetail; //Forward declaration

/**
 * @brief A mapper for linked reads from a PairedReadsDatastore.
 *
 * Supports partial remapping of unmapped reads or of a selection list.
 */

class PairedReadMapper {
public:
    PairedReadMapper(SequenceGraph &_sg, PairedReadsDatastore &_datastore) : sg(_sg),datastore(_datastore){
        reads_in_node.resize(sg.nodes.size());
    };
    void write(std::ofstream & output_file);
    void read(std::ifstream & input_file);
    void map_reads(std::unordered_set<uint64_t> const &  reads_to_remap={});
    void remap_all_reads();
    void map_read(uint64_t readID);
    void remove_obsolete_mappings();
    std::vector<uint64_t> size_distribution();
    void populate_orientation();
    /*void remap_reads();
    uint64_t process_reads_from_file(uint8_t, uint16_t, std::unordered_map<uint64_t , graphPosition> &, std::string , uint64_t, bool tags=false, std::unordered_set<uint64_t> const & reads_to_remap={});
    void save_to_disk(std::string filename);
    void load_from_disk(std::string filename);*/
    void print_stats();



    const SequenceGraph & sg;
    const PairedReadsDatastore & datastore;
    std::vector<std::vector<ReadMapping>> reads_in_node;
    std::vector<sgNodeID_t> read_to_node;//id of the main node if mapped, set to 0 to remap on next process
    //TODO: reading and writing this would simplify things??
    std::vector<bool> read_direction_in_node;//0-> fw, 1->rev;
    std::vector<uint64_t> rfdist;
    std::vector<uint64_t> frdist;
};

/**
 * @brief Analysis of all reads connecting two particular nodes.
 */

class PairedReadConnectivityDetail {
public:
    PairedReadConnectivityDetail(){};
    PairedReadConnectivityDetail(const PairedReadMapper & prm, sgNodeID_t source, sgNodeID_t dest);
    PairedReadConnectivityDetail& operator+=(const PairedReadConnectivityDetail& rhs){
        this->pairs_per_orientation[0] += rhs.pairs_per_orientation[0];
        this->pairs_per_orientation[1] += rhs.pairs_per_orientation[1];
        this->pairs_per_orientation[2] += rhs.pairs_per_orientation[2];
        this->pairs_per_orientation[3] += rhs.pairs_per_orientation[3];
        return *this;
    }

    uint64_t pairs_per_orientation[4]={0,0,0,0};
};

#endif //BSG_PAIREDREADMAPPER_HPP
