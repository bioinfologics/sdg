//
// Created by Ben Ward (EI) on 29/01/2018.
//

#ifndef BSG_SEQUENCEMAPPER_H
#define BSG_SEQUENCEMAPPER_H

#include "SequenceGraph.hpp"
#include <iostream>

typedef uint64_t seqID_t;

struct graphPosition{
    sgNodeID_t node;
    uint32_t pos;
};

enum MappingDirection {Nowhere, Forward, Backwards};

class SequenceMapping {
private:
    seqID_t seq_id;
    int32_t first_seq_pos;
    int32_t last_seq_pos;
    sgNodeID_t node;
    int32_t first_node_pos;
    int32_t last_node_pos;
    int32_t unique_matches;

    bool direction_will_continue(int32_t next_position) const;

public:
    SequenceMapping();
    bool operator==(const SequenceMapping &other);
    bool operator<(const SequenceMapping &other) const;
    friend std::ostream& operator<<(std::ostream& stream, const SequenceMapping& sm);
    void initiate_mapping(uint64_t sequence_id);
    bool ismatched();
    void start_new_mapping(sgNodeID_t nodeid, int32_t nodepos, int32_t seqpos);
    void extend(int32_t nodepos, int32_t seqpos);
    sgNodeID_t absnode() const;
    int32_t n_unique_matches() const;
    MappingDirection node_direction() const;
    MappingDirection seq_direction() const;
    bool mapping_continues(sgNodeID_t next_node, int32_t next_position) const;
};

class SequenceMapper {
    typedef std::unordered_map<uint64_t, graphPosition> UniqueKmerIndex;

private:
    SequenceGraph& sg;
    uint8_t k;
    UniqueKmerIndex graph_kmer_index;
    std::string query_seq_file;
    std::string output_prefix;
    uint64_t memory_limit;
    std::vector<std::vector<SequenceMapping>> mappings_in_node;
    std::unordered_map<seqID_t, std::vector<SequenceMapping>> mappings_of_sequence;

    void generate_unique_kmer_index(uint8_t k);
    void map_sequences_from_file(uint64_t min_matches, const std::string& filename);

public:
    SequenceMapper(SequenceGraph &_sg, uint8_t _k = 31) : sg(_sg), k(_k) {
        generate_unique_kmer_index(k);
        mappings_in_node.resize(sg.nodes.size());
    }

    void map_sequences(uint64_t min_matches, const std::string& filename) {
        map_sequences_from_file(min_matches, filename);
    }

    void print_mappings() const {
        auto v = mappings_of_sequence.find(0);
        if(v != mappings_of_sequence.end()){
            for (const auto &sm:v->second) {
                std::cout << sm << std::endl;
            }
        }
    }
};


#endif //BSG_SEQUENCEMAPPER_H
