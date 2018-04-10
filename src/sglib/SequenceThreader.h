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
#include <sglib/factories/KMerIDXFactory.h>
#include <sglib/logger/OutputLog.h>

typedef uint64_t seqID_t;

enum MappingDirection {Nowhere, Forward, Backwards};

class SequenceMapping;
class SequenceMappingThread;
class BridgedMappingThreads;

class SequenceThreader {
    using SequenceMappingStore = std::unordered_map<seqID_t, std::vector<SequenceMapping>>;
    using SequenceMappingPathsStore = std::unordered_map<seqID_t, std::vector<SequenceMappingThread>>;
    using BridgedMappingPathsStore = std::unordered_map<seqID_t, std::vector<BridgedMappingThreads>>;

public:
    explicit SequenceThreader(SequenceGraph &_sg, uint8_t _k = 31) : sg(_sg), k(_k), graph_kmer_index(sg, _k) {}
    // Kmer mapping.
    void map_sequences(const std::string &filename, const std::string &output) {
        output_prefix = output;
        query_seq_file = filename;

        map_sequences_from_file(filename);
    }
    void print_mappings(std::ostream& out, bool use_oldnames = false) const;

    // Connecting kmer mappings into threads through graph.
    void filter_mappings(double score = 90);
    void thread_mappings();
    void bridge_threads();
    void print_paths(std::ostream& out, bool use_oldnames = false) const;

    // File output and data dumping.
    void graph_threads_to_fasta(std::ofstream& output_file, bool use_oldnames = true) const;
    void query_threads_to_fasta(std::ofstream& output_file, bool use_oldnames = true) const;
    void bridged_graph_threads_to_fasta(std::ofstream& output_file, bool use_oldnames = true) const;
    void print_dark_nodes(std::ofstream& output_file) const;
    void print_full_node_diagnostics(std::ofstream& output_file) const;
    void print_unmapped_nodes(std::ofstream& output_file) const;

private:
    SequenceGraph& sg;
    UniqueKmerIndex graph_kmer_index;

    // Settings storage
    uint8_t k;
    std::string query_seq_file;
    std::string output_prefix;
    uint64_t memory_limit;

    // Results storage
    SequenceMappingStore mappings_of_sequence;
    SequenceMappingStore filtered_mappings_of_sequence;
    SequenceMappingPathsStore mapping_threads_of_sequence;
    BridgedMappingPathsStore  bridged_mapping_threads_of_sequence;
    std::vector<KmerIDX> unmapped_kmers;

    void map_sequences_from_file(const std::string &filename);
    std::vector<SequenceGraphPath> collect_paths(const sgNodeID_t seed, const sgNodeID_t target, unsigned int size_limit, unsigned int edge_limit);
};


class SequenceMapping {
    friend void SequenceThreader::print_mappings(std::ostream &out, bool use_oldnames) const;

public:
    SequenceMapping();
    bool operator==(const SequenceMapping &other);
    bool operator<(const SequenceMapping &other) const;
    friend std::ostream& operator<<(std::ostream& stream, const SequenceMapping& sm);
    void initiate_mapping(uint64_t sequence_id);
    bool ismatched();
    void start_new_mapping(const graphPosition& gpos, uint32_t seqpos, const UniqueKmerIndex& counts);
    void extend(int32_t nodepos, uint32_t seqpos);
    sgNodeID_t absnode() const;
    sgNodeID_t dirnode() const;
    int32_t n_unique_matches() const;
    MappingDirection node_direction() const;
    MappingDirection seq_direction() const;
    bool mapping_continues(const graphPosition& gpos) const;
    uint32_t query_start() const { return first_seq_pos; };
    uint32_t query_end() const { return last_seq_pos; };
    double match_score() const { return (double(matched_unique_kmers) / possible_unique_matches) * 100; };

private:
    seqID_t seq_id;
    uint32_t first_seq_pos;
    uint32_t last_seq_pos;
    sgNodeID_t node;
    int32_t first_node_pos;
    int32_t last_node_pos;
    int32_t matched_unique_kmers;
    uint64_t possible_unique_matches;
    uint64_t n_kmers_in_node;
    bool direction_will_continue(int32_t next_position) const;
};

class SequenceMappingThread {
public:

    explicit SequenceMappingThread(SequenceGraph& _sg) : ordered_mappings({}), node_path(_sg) {};
    SequenceMappingThread(const SequenceMappingThread& smt) = default;

    bool append_mapping(SequenceMapping mapping);
    void clear();

    size_t size() const;
    uint32_t query_start() const;
    uint32_t query_end() const;
    SequenceMapping first_mapping() const;
    SequenceMapping last_mapping() const;
    SequenceGraphPath get_graph_path() const { return node_path; }

    void print_path_header(std::ostream& output_file, bool use_oldnames = true) const;
    void print_sequence(std::ofstream& output_file) const;

private:
    std::vector<SequenceMapping> ordered_mappings;
    SequenceGraphPath node_path;
};

class BridgedMappingThreads {
public:
    BridgedMappingThreads(SequenceGraph& s, const SequenceMappingThread& smt)
            : sg(s), mapping_threads({}), bridging_paths({})
    {
        mapping_threads.emplace_back(smt);
    }

    BridgedMappingThreads& operator=(BridgedMappingThreads other);

    SequenceMapping last_mapping();

    uint32_t query_start() const;
    uint32_t query_end() const;

    void add_mapping_thread(const SequenceMappingThread& smt);

    int bridge_to_thread(const SequenceGraphPath& sgp, const SequenceMappingThread& smt);

    std::string get_complete_sequence() const;

    SequenceGraphPath get_complete_path() const;

private:
    SequenceGraph& sg;
    std::vector<SequenceMappingThread> mapping_threads;
    std::vector<SequenceGraphPath> bridging_paths;
};

#endif //BSG_SEQUENCEMAPPER_H
