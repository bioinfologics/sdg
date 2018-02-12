//
// Created by Ben Ward (EI) on 29/01/2018.
//

#ifndef BSG_SEQUENCEMAPPER_H
#define BSG_SEQUENCEMAPPER_H

#include "SequenceGraph.hpp"
#include <iostream>
#include <fstream>
#include <vector>
#include "UniqueKmerIndex.h"
#include "GraphDrawer.h"

typedef uint64_t seqID_t;

enum MappingDirection {Nowhere, Forward, Backwards};

class SequenceMapping {
private:
    seqID_t seq_id;
    int32_t first_seq_pos;
    int32_t last_seq_pos;
    sgNodeID_t node;
    int32_t first_node_pos;
    int32_t last_node_pos;
    int32_t matched_unique_kmers;
    uint64_t possible_unique_matches;
    uint64_t n_kmers_in_node;

    bool direction_will_continue(int32_t next_position) const;

public:
    SequenceMapping();
    bool operator==(const SequenceMapping &other);
    bool operator<(const SequenceMapping &other) const;
    friend std::ostream& operator<<(std::ostream& stream, const SequenceMapping& sm);
    void initiate_mapping(uint64_t sequence_id);
    bool ismatched();
    void start_new_mapping(const graphPosition& gpos, int32_t seqpos, const UniqueKmerIndex& counts);
    void extend(int32_t nodepos, int32_t seqpos);
    sgNodeID_t absnode() const;
    sgNodeID_t dirnode() const;
    int32_t n_unique_matches() const;
    MappingDirection node_direction() const;
    MappingDirection seq_direction() const;
    bool mapping_continues(const graphPosition& gpos) const;
    double POPUKM() const;
};

class SequenceMapper {
    typedef std::unordered_map<seqID_t, std::vector<SequenceMapping>> SequenceMappingStore;
    typedef std::unordered_map<seqID_t, std::vector<std::vector<SequenceMapping>>> SequenceMappingPathsStore;

private:
    SequenceGraph& sg;
    UniqueKmerIndex graph_kmer_index;

    uint8_t k;
    std::string query_seq_file;
    std::string output_prefix;
    uint64_t memory_limit;
    SequenceMappingStore mappings_of_sequence;
    SequenceMappingPathsStore paths_of_mappings_of_sequence;

    void map_sequences_from_file(uint64_t min_matches, const std::string& filename);

public:
    SequenceMapper(SequenceGraph &_sg, uint8_t _k = 31) : sg(_sg), k(_k), graph_kmer_index(sg, _k) {}

    void map_sequences(uint64_t min_matches, const std::string& filename, const std::string& output) {
        output_prefix = output;
        query_seq_file = filename;
        map_sequences_from_file(min_matches, filename);
    }

    void print_mappings() const {
        for( auto it = mappings_of_sequence.begin(); it != mappings_of_sequence.end(); ++it) {
            std::cout << "Mappings for query sequence: " << it->first << ":" << std::endl;
            for (const auto &sm:it->second) {
                std::cout << sm << std::endl;
            }
        }
    }

    void mappings_paths();

    void print_paths() const {
        for( auto it = paths_of_mappings_of_sequence.begin(); it != paths_of_mappings_of_sequence.end(); ++it) {
            std::cout << "Paths for query sequence: " << it->first << ":" << std::endl;
            for (const auto &path:it->second) {
                std::cout << "Path of " << path.size() << " mappings: ";
                for (const auto &sm:path) {
                    sgNodeID_t dirnode = sm.node_direction() == Forward ? sm.absnode() : -sm.absnode();
                    std::cout << dirnode << ", ";
                }
                std::cout << std::endl;
            }
        }
    }

    void paths_to_fasta(std::ofstream& output_file) const;
    void paint_paths(GraphDrawer& gd) const;
};


#endif //BSG_SEQUENCEMAPPER_H
