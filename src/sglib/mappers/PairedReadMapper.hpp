//
// Created by Bernardo Clavijo (EI) on 12/05/2018.
//

#ifndef BSG_PAIREDREADMAPPER_HPP
#define BSG_PAIREDREADMAPPER_HPP

#include <map>

#include "sglib/mappers/ReadMapping.hpp"
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
    void print_stats();

    std::vector<uint64_t> get_node_readpairs_ids(sgNodeID_t);

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
        this->orientation_paircount[0] += rhs.orientation_paircount[0];
        this->orientation_paircount[1] += rhs.orientation_paircount[1];
        this->orientation_paircount[2] += rhs.orientation_paircount[2];
        this->orientation_paircount[3] += rhs.orientation_paircount[3];
        return *this;
    }

    uint64_t orientation_paircount[4]={0,0,0,0};
};

#endif //BSG_PAIREDREADMAPPER_HPP
