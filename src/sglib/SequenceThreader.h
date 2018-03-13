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

class SequenceMapping;
class SequenceMappingThread;

enum MappingDirection {Nowhere, Forward, Backwards};

class MappingSummary {
public:
    bool operator<(MappingSummary rhs) {
        return id < rhs.id;
    }
private:
    unsigned long mapped_kmers, unmapped_kmers;
    seqID_t id;
    std::string name;
};

class SequenceThreader {
    typedef std::unordered_map<seqID_t, std::vector<SequenceMapping>> SequenceMappingStore;
    typedef std::unordered_map<seqID_t, std::vector<SequenceMappingThread>> SequenceMappingPathsStore;

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
    void thread_mappings();
    void connect_threads();
    void print_paths(std::ostream& out, bool use_oldnames = false) const;

    // File output and data dumping.
    void graph_threads_to_fasta(std::ofstream& output_file, bool use_oldnames = true) const;
    void query_threads_to_fasta(std::ofstream& output_file, bool use_oldnames = true) const;
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
    SequenceMappingPathsStore mapping_threads_of_sequence;
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
    void start_new_mapping(const graphPosition& gpos, int32_t seqpos, const UniqueKmerIndex& counts);
    void extend(int32_t nodepos, int32_t seqpos);
    sgNodeID_t absnode() const;
    sgNodeID_t dirnode() const;
    int32_t n_unique_matches() const;
    MappingDirection node_direction() const;
    MappingDirection seq_direction() const;
    bool mapping_continues(const graphPosition& gpos) const;
    uint32_t query_start() const {return first_seq_pos;};
    uint32_t query_end() const {return last_seq_pos;};

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

    bool append_mapping(SequenceMapping mapping) {
        auto dn = mapping.dirnode();
        sglib::OutputLog(sglib::LogLevels::DEBUG) << "Trying to add to path as: " << dn << std::endl;
        auto path_success = node_path.append_to_path(dn);
        if (path_success) {
            sglib::OutputLog(sglib::LogLevels::DEBUG) << "Was able to append " << dn << " to the path." << std::endl;
            ordered_mappings.emplace_back(mapping);
        } else {
            sglib::OutputLog(sglib::LogLevels::DEBUG) << "Was not able to append " << dn << " to the path." << std::endl;
        }
        return path_success;
    }

    void clear() {
        ordered_mappings.clear();
        node_path.clear();
    }

    size_t size() const {
        return ordered_mappings.size();
    }

    uint32_t query_start() const {
        return ordered_mappings.front().query_start();
    }

    uint32_t query_end() const {
        return ordered_mappings.back().query_end();
    }

    void print_path_header(std::ofstream& output_file, bool use_oldnames = true) const {
        auto fasta_header = node_path.get_fasta_header(use_oldnames);
        fasta_header.erase(fasta_header.begin());
        output_file << fasta_header;
    }

    void print_sequence(std::ofstream& output_file) const {
        auto s = node_path.get_sequence();
        output_file << s;
    }

    SequenceMapping first_mapping() const {
        return ordered_mappings.front();
    }

    SequenceMapping last_mapping() const {
        return ordered_mappings.back();
    }

private:
    std::vector<SequenceMapping> ordered_mappings;
    SequenceGraphPath node_path;
};




#endif //BSG_SEQUENCEMAPPER_H
