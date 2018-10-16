//
// Created by Bernardo Clavijo (EI) on 11/02/2018.
//

#ifndef BSG_LINKEDREADMAPPER_HPP
#define BSG_LINKEDREADMAPPER_HPP

#include <map>

#include "sglib/graph/SequenceGraph.hpp"
#include <sglib/datastores/LinkedReadsDatastore.hpp>
#include <sglib/types/MappingTypes.hpp>
#include <sglib/indexers/UniqueKmerIndex.hpp>

class UniqueKmerIndex;
class Unique63merIndex;
class TagNeighbour {
public:
    TagNeighbour(){};
    TagNeighbour(sgNodeID_t n, float s):node(n),score(s){};
    sgNodeID_t node;
    float score; //breaking the latte principle
};

/**
 * @brief A mapper for linked reads from a LinkedReadsDatastore.
 *
 * Supports partial remapping of unmapped reads or of a selection list.
 */
class LinkedReadMapper {

public:
    LinkedReadMapper(SequenceGraph &_sg, LinkedReadsDatastore &_datastore, const UniqueKmerIndex &uki, const Unique63merIndex &u63i) :
    sg(_sg),
    datastore(_datastore),
    kmer_to_graphposition(uki),
    k63mer_to_graphposition(u63i)
    {
        reads_in_node.resize(sg.nodes.size());
    };
    void write(std::ofstream & output_file);
    void read(std::ifstream & input_file);
    void map_reads(std::unordered_set<uint64_t> const &  reads_to_remap={});
    void remap_all_reads();
    LinkedReadMapper operator=(const LinkedReadMapper &other);


    void map_reads63(std::unordered_set<uint64_t> const &  reads_to_remap={});
    void remap_all_reads63();

    void remove_obsolete_mappings();
    /*void remap_reads();
    uint64_t process_reads_from_file(uint8_t, uint16_t, std::unordered_map<uint64_t , graphPosition> &, std::string , uint64_t, bool tags=false, std::unordered_set<uint64_t> const & reads_to_remap={});
    void save_to_disk(std::string filename);
    void load_from_disk(std::string filename);*/
    void print_stats(){};
    std::set<bsg10xTag> get_node_tags(sgNodeID_t n);
    std::map<bsg10xTag, std::vector<sgNodeID_t>> get_tag_nodes(uint32_t min_nodes = 2,
                                                               const std::vector<bool> &selected_nodes = {});
    std::vector<std::pair<sgNodeID_t , sgNodeID_t >> get_tag_neighbour_nodes(uint32_t min_shared,const std::vector<bool> & selected_nodes={});
    void compute_all_tag_neighbours(int min_size,float min_score);
    void write_tag_neighbours(std::string filename);
    void read_tag_neighbours(std::string filename);

    SequenceGraph & sg;
    const UniqueKmerIndex& kmer_to_graphposition;
    const Unique63merIndex& k63mer_to_graphposition;
    LinkedReadsDatastore &datastore;
    std::vector<std::vector<ReadMapping>> reads_in_node;
    std::vector<sgNodeID_t> read_to_node;//id of the main node if mapped, set to 0 to remap on next process
    std::vector<std::vector<TagNeighbour>> tag_neighbours; //not persisted yet!

    static const bsgVersion_t min_compat;
};
#endif //BSG_LINKEDREADMAPPER_HPP
