//
// Created by Ben Ward (EI) on 29/01/2018.
//

#ifndef BSG_SEQUENCEMAPPER_H
#define BSG_SEQUENCEMAPPER_H

#include <iostream>
#include <fstream>
#include <vector>
#include <set>
#include <sglib/SequenceGraph.h>
#include <sglib/UniqueKmerIndex.h>

typedef uint64_t seqID_t;

class SequenceMapping;

enum MappingDirection {Nowhere, Forward, Backwards};

class SequenceThreader {
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
    std::set<SequenceGraphPath> all_unique_paths() const;
    void print_mapping_path_name(const std::vector<SequenceMapping>& path, std::ofstream& output_file) const;

public:
    explicit SequenceThreader(SequenceGraph &_sg, uint8_t _k = 31) : sg(_sg), k(_k), graph_kmer_index(sg, _k) {}

    void map_sequences(uint64_t min_matches, const std::string& filename, const std::string& output) {
        output_prefix = output;
        query_seq_file = filename;

        map_sequences_from_file(min_matches, filename);
    }

    void print_mappings(std::ostream& out, bool use_oldnames = false) const;
    void mappings_paths();
    void print_paths(std::ostream& out, bool use_oldnames = false) const;
    void graph_paths_to_fasta(std::ofstream& output_file) const;
    void query_paths_to_fasta(std::ofstream& output_file) const;
    void print_unique_paths_sizes(std::ofstream& output_file) const;
    void print_dark_nodes(std::ofstream& output_file) const;
};

class SequenceMapping {
    friend void SequenceThreader::query_paths_to_fasta(std::ofstream& output_file) const;
    friend void SequenceThreader::print_mappings(std::ostream &out, bool use_oldnames) const;
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

#endif //BSG_SEQUENCEMAPPER_H
