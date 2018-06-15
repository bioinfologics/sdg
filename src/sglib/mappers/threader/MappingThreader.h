//
// Created by Ben Ward (EI) on 21/05/2018.
//

#ifndef BSG_MAPPINGTHREADER_H
#define BSG_MAPPINGTHREADER_H

#include <sglib/mappers/threader/NodeMapper.h>

class MappingThread;

// This thing constructs and stored MappingThread objects from the NodeMappings constructed by and stored in a
// Node mapper.

class MappingThreader {

public:
    using SequenceMappingPathsStore = std::unordered_map<seqID_t, std::vector<MappingThread>>;
    MappingThreader(NodeMapper& mppr) : mapper(mppr) {};
    const std::string getQueryFile() const { return mapper.getQueryFile(); }
    const SequenceGraph& getGraph() const { return mapper.getGraph(); }
    const uint8_t getK() const { return mapper.getK(); }
    const SequenceMappingPathsStore& getThreads() const { return mapping_threads_of_sequence; }
    void thread_mappings();
    void graph_threads_to_fasta(std::ofstream& output_file, bool use_oldnames) const;
    void query_threads_to_fasta(std::ofstream& output_file, bool use_oldnames) const;
    void calculate_reference_inclusion();

private:
    NodeMapper &mapper;
    SequenceMappingPathsStore mapping_threads_of_sequence;

    bool append_mapping_trivially(MappingThread &thread, NodeMapping mapping);
    bool append_mapping_with_path(MappingThread &thread, NodeMapping mapping, SequenceGraphPath &path);
};


class MappingThread {
    friend class MappingThreader;

public:
    explicit MappingThread(const SequenceGraph& _sg) : ordered_mappings(), node_path(_sg) {};
    MappingThread(const MappingThread& smt) = default;
    void clear() { ordered_mappings.clear(); node_path.clear(); };
    size_t size() const { return ordered_mappings.size(); };
    uint32_t queryStart() const { return firstMapping().query_start(); }
    uint32_t queryEnd() const { return lastMapping().query_end(); }
    uint32_t querySize() const;
    NodeMapping firstMapping() const { return ordered_mappings.front(); }
    NodeMapping lastMapping() const { return ordered_mappings.back(); }
    SequenceGraphPath get_graph_path() const { return node_path; }
    void print_path_header(std::ostream& output_file, bool use_oldnames = true) const;
    std::string get_sequence() const { return node_path.get_sequence(); };
    void print_sequence(std::ofstream& output_file) const { output_file << get_sequence(); };


private:
    std::vector<NodeMapping> ordered_mappings;
    SequenceGraphPath node_path;

    // Private methods for conveinience of MappingThreader friend.
    std::vector<NodeMapping>& orderedMappings() { return ordered_mappings; }
    SequenceGraphPath& nodePath() { return node_path; }
};

#endif //BSG_MAPPINGTHREADER_H
