//
// Created by Bernardo Clavijo (EI) on 11/02/2018.
//

#ifndef BSG_LINKEDREADMAPPER_HPP
#define BSG_LINKEDREADMAPPER_HPP

#include <map>

#include "sglib/SequenceGraph.h"
#include <sglib/datastores/LinkedReadsDatastore.hpp>


/**
 * @brief A mapper for linked reads from a LinkedReadsDatastore.
 *
 * Supports partial remapping of unmapped reads or of a selection list.
 */
class LinkedReadMapper {
public:
    LinkedReadMapper(SequenceGraph &_sg, LinkedReadsDatastore &_datastore) : sg(_sg),datastore(_datastore){
        reads_in_node.resize(sg.nodes.size());
    };
    void update_graph_index();
    void map_reads(std::unordered_set<uint64_t> const &  reads_to_remap={});
    void map_read(uint64_t readID);
    void remove_obsolete_mappings();
    /*void remap_reads();
    uint64_t process_reads_from_file(uint8_t, uint16_t, std::unordered_map<uint64_t , graphPosition> &, std::string , uint64_t, bool tags=false, std::unordered_set<uint64_t> const & reads_to_remap={});
    void save_to_disk(std::string filename);
    void load_from_disk(std::string filename);*/
    void print_stats(){};

    SequenceGraph & sg;
    LinkedReadsDatastore &datastore;
    std::unordered_map<uint64_t, graphPosition> kmer_to_graphposition;
    uint64_t memlimit;
    std::vector<std::vector<ReadMapping>> reads_in_node;
    std::vector<sgNodeID_t> read_to_node;//id of the main node if mapped, set to 0 to remap on next process
};


#endif //BSG_LINKEDREADMAPPER_HPP
