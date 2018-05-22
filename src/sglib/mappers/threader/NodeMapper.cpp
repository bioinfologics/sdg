//
// Created by Ben Ward (EI) on 18/05/2018.
//

#include <sglib/mappers/threader/NodeMapper.h>

// Public Methods...

NodeMapper::NodeMapper(SequenceGraph &sg, uniqueKmerIndex &uki) : sg(sg), graph_kmer_index(uki) {}

void NodeMapper::mapSequences(const std::string &filename) {
    query_seq_file = filename;
    map_sequences_from_file(filename);
}

void NodeMapper::write_mappings_to_binary_file(std::string filename) const {
    sglib::OutputLog() << "Dumping node mappings" << std::endl;
    std::ofstream outf(filename, std::ios_base::binary);
    auto mapSize(mappings_of_sequence.size());
    outf.write(reinterpret_cast<const char *>(&mapSize), sizeof(mapSize));
    for(const auto &it : mappings_of_sequence) {
        // Get the seqID of the current thing.
        auto sequence_id = it.first;
        outf.write(reinterpret_cast<const char *>(&sequence_id), sizeof(sequence_id));
        auto vecSize(it.second.size());
        outf.write(reinterpret_cast<const char *>(&vecSize), sizeof(vecSize));
        outf.write(reinterpret_cast<const char *>(it.second.data()), vecSize * sizeof(NodeMapping));
    }
}

void NodeMapper::read_mappings_from_binary_file(std::string filename) {
    sglib::OutputLog() << "Loading node mappings" << std::endl;
    std::ifstream inf(filename, std::ios_base::binary);
    auto mapSize(mappings_of_sequence.size());
    inf.read(reinterpret_cast<char *>(&mapSize), sizeof(mapSize));
    mappings_of_sequence.reserve(mapSize);
    seqID_t seqid;
    std::vector<NodeMapping> vec;
    auto vecSize { vec.size() };
    for (decltype(mapSize) n = 0; n < mapSize; ++n){
        inf.read(reinterpret_cast<char *>(&seqid), sizeof(seqid));
        inf.read(reinterpret_cast<char *>(&vecSize), sizeof(vecSize));
        vec.resize(vecSize);
        inf.read(reinterpret_cast<char *>(vec.data()), vecSize * sizeof(NodeMapping));
        mappings_of_sequence.emplace(seqid, vec);
    }
}

void NodeMapper::filterMappings(double score) {
    sglib::OutputLog(sglib::LogLevels::INFO) << "Filtering mappings using a threshold of " << score << "%" << std::endl;
    for (auto& sequence_mappings : mappings_of_sequence) {
        auto start { sequence_mappings.second.begin() };
        auto end { sequence_mappings.second.end() };
        auto newend {std::remove_if(start, end, [score](NodeMapping sm) -> bool { return sm.match_score() < score; })};
        sequence_mappings.second.erase(newend, sequence_mappings.second.end());
    }
}

void NodeMapper::print_unmapped_nodes(std::ofstream& output_file) const {
    output_file << "nodeID\tseq_size\tn_kmers\tn_unique_kmers\tn_forward_links\tn_backward_links" << std::endl;
    std::set<sgNodeID_t> mapped, all, diff;
    for (sgNodeID_t i = 0; i < sg.nodes.size(); i++) {
        all.insert(i);
    }
    for (const auto& sm : mappings_of_sequence) {
        for (const auto& m : sm.second) {
            mapped.insert(m.absnode());
        }
    }
    std::set_difference(all.begin(), all.end(), mapped.begin(), mapped.end(), std::inserter(diff, diff.end()));
    for (const auto& node : diff) {
        output_file << sg.nodeID_to_name(node) << '\t' << sg.nodes[node].sequence.size() << '\t';
        output_file << std::to_string(graph_kmer_index.total_kmers_in_node(node)) << '\t';
        output_file << std::to_string(graph_kmer_index.unique_kmers_in_node(node)) << '\t';
        output_file << std::to_string(sg.get_fw_links(node).size()) << '\t';
        output_file << std::to_string(sg.get_bw_links(node).size()) << '\t';
        output_file << std::endl;
    }
}
/*
void NodeMapper::printMappings(std::ofstream& outputFile) const {
    for (const auto& sm : mappings_of_sequence) {
        for (const auto& m : sm.second) {
            outputFile << m << std::endl;
        }
    }
}
*/
// Private Methods...

void NodeMapper::start_new_mapping(NodeMapping &mapping, const graphStrandPos gpos, uint32_t seqpos) {
    mapping.first_seq_pos = seqpos;
    mapping.last_seq_pos = seqpos;
    mapping.node = gpos.node;
    mapping.first_node_pos = gpos.pos;
    mapping.last_node_pos = gpos.pos;
    mapping.matched_unique_kmers = 1;
    mapping.possible_unique_matches = graph_kmer_index.unique_kmers_in_node(gpos.node);
    mapping.n_kmers_in_node = graph_kmer_index.total_kmers_in_node(gpos.node);
}

void NodeMapper::map_sequences_from_file(const std::string &filename) {
    sglib::OutputLog(sglib::LogLevels::INFO) << "Mapping sequence kmers to graph" << std::endl;
    FastaReader<FastaRecord> fastaReader({0}, filename);
    std::atomic<uint64_t> sequence_count(0);

#pragma omp parallel shared(fastaReader, mappings_of_sequence, sequence_unmapped_kmers, query_seq_sizes, query_seq_names)
    {
        FastaRecord sequence;
        bool c;

#pragma omp critical (getrecord)
        {
            c = fastaReader.next_record(sequence);
        }

        while (c) {
            sequence_count++;
            std::vector<NodeMapping> mappings;
            std::vector<KmerIDX> unmapped_kmers;
            sglib::OutputLog(sglib::LogLevels::INFO) << "Mapping sequence: " << sequence.name << " (" << sequence.id << ')' << std::endl;
            std::tie(mappings, unmapped_kmers) = map_sequence_to_graph(sequence);

#pragma omp critical (store)
            {
                sglib::OutputLog(sglib::LogLevels::INFO) << "Failed to map " << unmapped_kmers.size() << " kmers from the reference" << std::endl;
                mappings_of_sequence[sequence.id] = mappings;
                sequence_unmapped_kmers[sequence.id] = unmapped_kmers;
                query_seq_sizes[sequence.id] = sequence.seq.size();
                query_seq_names[sequence.id] = sequence.name;
            }

#pragma omp critical (getrecord)
            {
                c = fastaReader.next_record(sequence);
            }
        }
    }
}

std::tuple<std::vector<NodeMapping>, std::vector<KmerIDX>> NodeMapper::map_kmers_to_graph(seqID_t id, std::vector<KmerIDX>& kmers) {
    NodeMapping mapping;
    mapping.initiate_mapping(id);
    // For each kmer on sequence.
    uint64_t mapped_kmers_count {0};
    std::vector<NodeMapping> sequence_mappings;
    std::vector<KmerIDX> unmapped_kmers;
    for (auto &sk:kmers) {
        bool found_kmer;
        graphStrandPos graph_pos;
        std::tie(found_kmer, graph_pos) = graph_kmer_index.find_unique_kmer_in_graph(sk.kmer);
        // IF KMER EXISTS ON GRAPH
        if (found_kmer) {
            sglib::OutputLog(sglib::LogLevels::DEBUG) << "Found kmer: " << sk.kmer << "in graph" << std::endl;
            sglib::OutputLog(sglib::LogLevels::DEBUG) << '(' << graph_pos.node << ", " << graph_pos.pos << ')' << std::endl;
            mapped_kmers_count++;
            // IF THE KMER MATCH IS THE FIRST MATCH FOR THE MAPPING...
            if (!mapping.ismatched()) {
                sglib::OutputLog(sglib::LogLevels::DEBUG) << "Kmer match is the first for the mapping." << std::endl;
                start_new_mapping(mapping, graph_pos, sk.pos);
            } // IF THE KMER MATCH IS NOT THE FIRST MATCH FOR THE MAPPING...
            else {
                sglib::OutputLog(sglib::LogLevels::DEBUG) << "Kmer match is not the first match for the mapping." << std::endl;
                // THERE ARE TWO SITUATIONS WHERE WE WOULD START A NEW MAPPING...
                // A). WE ARE MAPPING TO A COMPLETELY DIFFERENT NODE...
                // B). THE DIRECTION OF THE MAPPING IS FLIPPED...

                if (!mapping.mapping_continues(graph_pos)) {
                    sglib::OutputLog(sglib::LogLevels::DEBUG) << "Kmer match does not continues the current mapping." << std::endl;
                    sequence_mappings.push_back(mapping);
                    start_new_mapping(mapping, graph_pos, sk.pos);
                } else { // IF NODE KMER MAPS TO IS THE SAME, EXTEND THE CURRENT MAPPING...
                    sglib::OutputLog(sglib::LogLevels::DEBUG) << "Kmer match continues the current mapping." << std::endl;
                    mapping.extend(graph_pos.pos, sk.pos);
                }
            }
        } else {
            unmapped_kmers.emplace_back(sk.kmer);
        }
    }
    // TODO : Check if this last push_back is required
    if (mapping.ismatched()) sequence_mappings.push_back(mapping);
    sglib::OutputLog(sglib::LogLevels::INFO) << "Mapped " << mapped_kmers_count << " kmers from the reference." << std::endl;
    return std::make_tuple(std::move(sequence_mappings), std::move(unmapped_kmers));
}

std::tuple<std::vector<NodeMapping>, std::vector<KmerIDX>> NodeMapper::map_sequence_to_graph(FastaRecord& seq) {
    uint8_t k = graph_kmer_index.get_k();
    kmerIDXFactory<FastaRecord> kf({k});
    std::vector<KmerIDX> sequence_kmers(0);
    kf.setFileRecord(seq);
    kf.next_element(sequence_kmers);
    return map_kmers_to_graph(seq.id, sequence_kmers);
}

// NodeMapping...

// Public Methods...

NodeMapping::NodeMapping(){
    // Just clean the structure, OSX doesn't give you clean memory
    memset(this, 0, sizeof(NodeMapping));
}

NodeMapping::NodeMapping(seqID_t id, uint32_t seq_first, uint32_t seq_last, sgNodeID_t node, int32_t node_first,
                         int32_t node_last, int32_t muk, uint64_t puk, uint64_t nk)
        : seq_id(id), first_seq_pos(seq_first), last_seq_pos(seq_last), node(node), first_node_pos(node_first),
          last_node_pos(node_last), matched_unique_kmers(muk), possible_unique_matches(puk), n_kmers_in_node(nk)
{}

bool NodeMapping::operator==(const NodeMapping &other) const {
    if (seq_id == other.seq_id && first_seq_pos == other.first_seq_pos && last_seq_pos == other.last_seq_pos &&
        node == other.node && first_node_pos == other.first_node_pos && last_node_pos == other.last_node_pos &&
        matched_unique_kmers == other.matched_unique_kmers && possible_unique_matches == other.possible_unique_matches &&
        n_kmers_in_node == other.n_kmers_in_node) {
        return true;
    } else {
        return false;
    }
};

bool NodeMapping::operator<(const NodeMapping &other) const {
    if (node != other.node) return node < other.node;
    return seq_id < other.seq_id;
};

void NodeMapping::initiate_mapping(seqID_t sequence_id) {
    seq_id = sequence_id;
    first_seq_pos = 0;
    last_seq_pos = 0;
    node = 0;
    first_node_pos = 0;
    last_node_pos = 0;
    matched_unique_kmers = 0;
    possible_unique_matches = 0;
}

bool NodeMapping::ismatched(){
    return node != 0;
}

void NodeMapping::extend(int32_t nodepos, uint32_t seqpos) {
    matched_unique_kmers++;
    last_node_pos = nodepos;
    last_seq_pos = seqpos;
}

sgNodeID_t NodeMapping::absnode() const {
    return std::abs(node);
}

sgNodeID_t NodeMapping::dirnode() const {
    return node_direction() == Forward ? absnode() : -absnode();
}

int32_t NodeMapping::n_unique_matches() const {
    return matched_unique_kmers;
}

MappingDirection NodeMapping::node_direction() const {
    auto d = last_node_pos - first_node_pos;
    return d > 0 ? Forward : d < 0 ? Backwards : Nowhere;
}

MappingDirection NodeMapping::seq_direction() const {
    auto d = last_seq_pos - first_seq_pos;
    return d > 0 ? Forward : d < 0 ? Backwards : Nowhere;
}

uint32_t NodeMapping::query_start() const {
    return first_seq_pos;
};

uint32_t NodeMapping::query_end() const {
    return last_seq_pos;
};

double NodeMapping::match_score() const {
    return (double(matched_unique_kmers) / possible_unique_matches) * 100;
};

std::ostream& operator<<(std::ostream& stream, const NodeMapping& sm) {
    auto sd = sm.seq_direction();
    auto nd = sm.node_direction();
    auto seqdir = sd == Forward ? "Forward" : sd == Backwards ? "Backwards" : "DIRERR";
    auto nodedir = nd == Forward ? "Forward" : nd == Backwards ? "Backwards" : "DIRERR";
    stream << "Sequence: " << sm.seq_id << " from: ";
    stream << sm.first_seq_pos << " : " << sm.last_seq_pos << " (" << seqdir << ")";
    stream << ", maps to node: " << sm.absnode();
    stream << " from: " << sm.first_node_pos << " : " << sm.last_node_pos  << " (" << nodedir << "). ";
    stream << sm.matched_unique_kmers << " / " << sm.possible_unique_matches << " unique kmers matched.";
    stream << " (" << sm.n_kmers_in_node << " kmers exist in node).";
    return stream;
}

// Private methods...

bool NodeMapping::direction_will_continue(int32_t next_position) const {
    auto direction = node_direction();
    if(direction == Forward) {
        return next_position > last_node_pos;
    } else if(direction == Backwards) {
        return next_position < last_node_pos;
    } else if(direction == Nowhere) {
        return true;
    }
}

bool NodeMapping::mapping_continues(const graphStrandPos& gpos) const {
    auto same_node = absnode() == std::abs(gpos.node);
    auto direction_continues = direction_will_continue(gpos.pos);
    return same_node and direction_continues;
}

