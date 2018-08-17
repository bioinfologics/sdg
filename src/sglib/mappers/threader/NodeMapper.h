//
// Created by Ben Ward (EI) on 18/05/2018.
//

#ifndef BSG_NODEMAPPER_H
#define BSG_NODEMAPPER_H

#include <cstdint>
#include <sglib/types/GenericTypes.hpp>
#include <sglib/types/KmerTypes.hpp>
#include <sglib/graph/SequenceGraph.hpp>
#include <sglib/indexers/UniqueKmerIndex.hpp>

class NodeMapper;

enum MappingDirection {Nowhere, Forward, Backwards};

class NodeMapping {
    friend NodeMapper;
    friend std::ostream& operator<<(std::ostream& stream, const NodeMapping& sm);

public:
    NodeMapping();
    NodeMapping(seqID_t id, uint32_t seq_first, uint32_t seq_last, sgNodeID_t node, int32_t node_first,
                int32_t node_last, int32_t muk, uint64_t puk, uint64_t nk);
    bool operator==(const NodeMapping &other) const;
    bool operator<(const NodeMapping &other) const;
    void initiate_mapping(seqID_t sequence_id);
    bool ismatched();
    //void start_new_mapping(const graphStrandPos& gpos, uint32_t seqpos, const uniqueKmerIndex& counts);
    void extend(int32_t nodepos, uint32_t seqpos);
    sgNodeID_t absnode() const;
    sgNodeID_t dirnode() const;
    int32_t n_unique_matches() const;
    MappingDirection node_direction() const;
    MappingDirection seq_direction() const;
    bool mapping_continues(const graphStrandPos& gpos) const;
    uint32_t query_start() const;
    uint32_t query_end() const;
    double match_score() const;

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

class NodeMapper {
    using SequenceUnmappedKmers = std::unordered_map<seqID_t, std::vector<KmerIDX>>;

public:
    using SequenceMappingStore = std::unordered_map<seqID_t, std::vector<NodeMapping>>;
    NodeMapper(SequenceGraph &sg, UniqueKmerIndex &uki);
    void mapSequences(const std::string &filename);
    void write_mappings_to_binary_file(std::string filename) const;
    void read_mappings_from_binary_file(std::string filename);
    void filterMappings(double score);
    void print_unmapped_nodes(std::ofstream& output_file) const;
    //void printMappings(std::ofstream& outputFile) const;

    void printMappings(std::ostream& out, bool use_oldnames) const {
        out << "ref_sequence\tref_start\tref_end\tgraph_node\tnode_start\tnode_end\tunique_kmers_matched\tn_unique_kmers\tn_kmers" << std::endl;
        for(const auto &it:mappings_of_sequence) {
            for (const auto &sm:it.second) {
                out << sm.seq_id << '\t' << sm.first_seq_pos << '\t' << sm.last_seq_pos << '\t';
                out << (use_oldnames ? sg.nodeID_to_name(sm.absnode()) : std::to_string(sm.absnode())) << '\t';
                out << sm.first_node_pos << '\t' << sm.last_node_pos  << '\t' << sm.matched_unique_kmers << '\t' << sm.possible_unique_matches;
                out << '\t' << sm.n_kmers_in_node << std::endl;
            }
        }
    }

    // Accessors - only return internal private variables that may be large by CONST ref.
    const SequenceMappingStore& getMappings() const { return mappings_of_sequence; }
    const SequenceGraph& getGraph() const { return sg; }
    const std::string getQueryFile() const { return query_seq_file; }
    const uint8_t getK() const { return graph_kmer_index.get_k(); }
    const std::unordered_map<seqID_t, std::string>& querySeqNames() { return query_seq_names; }
    const std::unordered_map<seqID_t, size_t>& querySeqSizes() { return query_seq_sizes; }

private:
    // Graph storage
    const SequenceGraph& sg;
    UniqueKmerIndex& graph_kmer_index;

    // Settings storage
    std::string query_seq_file;

    // Seq info storage
    std::unordered_map<seqID_t, size_t> query_seq_sizes;
    std::unordered_map<seqID_t, std::string> query_seq_names;

    // Mapping storage
    SequenceMappingStore mappings_of_sequence;
    SequenceUnmappedKmers sequence_unmapped_kmers;


    void start_new_mapping(NodeMapping &mapping, graphStrandPos gpos, uint32_t seqpos);
    void map_sequences_from_file(const std::string &filename);
    std::tuple<std::vector<NodeMapping>, std::vector<KmerIDX>> map_kmers_to_graph(seqID_t id, std::vector<KmerIDX>& kmers);
    std::tuple<std::vector<NodeMapping>, std::vector<KmerIDX>> map_sequence_to_graph(FastaRecord& seq);
};



#endif //BSG_NODEMAPPER_H
