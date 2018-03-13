//
// Created by Ben Ward (EI) on 29/01/2018.
//

#include <numeric>
#include <atomic>
#include <algorithm>
#include <stack>
#include <sglib/SequenceThreader.h>
#include <sglib/readers/FileReader.h>
#include <sglib/logger/OutputLog.h>

SequenceMapping::SequenceMapping(){
    // Just clean the structure, OSX doesn't give you clean memory
    bzero(this, sizeof(SequenceMapping));
}

bool SequenceMapping::operator==(const SequenceMapping &other){
    return this == &other;
};

bool SequenceMapping::operator<(const SequenceMapping &other) const {
    if (node != other.node) return node < other.node;
    return seq_id < other.seq_id;
};

std::ostream& operator<<(std::ostream& stream, const SequenceMapping& sm) {
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

void SequenceMapping::initiate_mapping(uint64_t sequence_id) {
    seq_id = sequence_id;
    first_seq_pos = 0;
    last_seq_pos = 0;
    node = 0;
    first_node_pos = 0;
    last_node_pos = 0;
    matched_unique_kmers = 0;
    possible_unique_matches = 0;
}

bool SequenceMapping::ismatched(){
    return node != 0;
}

void SequenceMapping::start_new_mapping(const graphPosition& gpos, int32_t seqpos, const UniqueKmerIndex& index) {
    first_seq_pos = seqpos;
    last_seq_pos = seqpos;
    node = gpos.node;
    first_node_pos = gpos.pos;
    last_node_pos = gpos.pos;
    matched_unique_kmers = 1;
    possible_unique_matches = index.unique_kmers_in_node(gpos.node);
    n_kmers_in_node = index.total_kmers_in_node(gpos.node);
}

void SequenceMapping::extend(int32_t nodepos, int32_t seqpos) {
    matched_unique_kmers++;
    last_node_pos = nodepos;
    last_seq_pos = seqpos;
}

sgNodeID_t SequenceMapping::absnode() const {
    return std::abs(node);
}

sgNodeID_t SequenceMapping::dirnode() const {
    return node_direction() == Forward ? absnode() : -absnode();
}

int32_t SequenceMapping::n_unique_matches() const {
    return matched_unique_kmers;
}

MappingDirection SequenceMapping::node_direction() const {
    auto d = last_node_pos - first_node_pos;
    return d > 0 ? Forward : d < 0 ? Backwards : Nowhere;
}

MappingDirection SequenceMapping::seq_direction() const {
    auto d = last_seq_pos - first_seq_pos;
    return d > 0 ? Forward : d < 0 ? Backwards : Nowhere;
}

bool SequenceMapping::direction_will_continue(int32_t next_position) const {
    auto direction = node_direction();
    //std::cout << "Node direction: " << direction << std::endl;
    if(direction == Forward) {
        return next_position > last_node_pos;
    } else if(direction == Backwards) {
        return next_position < last_node_pos;
    } else if(direction == Nowhere) {
        return true;
    }
}

bool SequenceMapping::mapping_continues(const graphPosition& gpos) const {
    auto same_node = absnode() == std::abs(gpos.node);
    auto direction_continues = direction_will_continue(gpos.pos);
    return same_node and direction_continues;
}

void SequenceThreader::map_sequences_from_file(const std::string &filename) {

    sglib::OutputLog(sglib::LogLevels::INFO) << "Mapping sequence kmers to graph." << std::endl;
    FastaReader<FastaRecord> fastaReader({0}, filename);
    std::atomic<uint64_t> mapped_kmers_count(0), sequence_count(0);

#pragma omp parallel shared(fastaReader, mappings_of_sequence, unmapped_kmers)
    {
        FastaRecord sequence;
        std::vector<KmerIDX> seqkmers;
        kmerIDXFactory<FastaRecord> kf({k});
        SequenceMapping mapping;
        bool c;

#pragma omp critical (readrec)
        {
            c = fastaReader.next_record(sequence);
        }

        while (c) {
            sequence_count++;
            mapping.initiate_mapping((uint64_t) sequence.id);

            // get all Kmers from sequence.
            seqkmers.clear();

            kf.setFileRecord(sequence);
            kf.next_element(seqkmers);

            // For each unique kmer on sequence.
            std::vector<SequenceMapping> sequence_mappings;
            for (auto &sk:seqkmers) {
                bool found_kmer;
                graphPosition graph_pos;
                std::tie(found_kmer, graph_pos) = graph_kmer_index.find_unique_kmer_in_graph(sk.kmer);
                // IF KMER EXISTS ON GRAPH
                if (found_kmer) {
                    sglib::OutputLog(sglib::LogLevels::DEBUG) << "Found kmer: " << sk.kmer << "in graph" << std::endl;
                    sglib::OutputLog(sglib::LogLevels::DEBUG) << '(' << graph_pos.node << ", " << graph_pos.pos << ')' << std::endl;
                    mapped_kmers_count++;
                    // IF THE KMER MATCH IS THE FIRST MATCH FOR THE MAPPING...
                    if (!mapping.ismatched()) {
                        sglib::OutputLog(sglib::LogLevels::DEBUG) << "Kmer match is the first for the mapping." << std::endl;
                        mapping.start_new_mapping(graph_pos, sk.pos, graph_kmer_index);
                    } // IF THE KMER MATCH IS NOT THE FIRST MATCH FOR THE MAPPING...
                    else {
                        sglib::OutputLog(sglib::LogLevels::DEBUG) << "Kmer match is not the first match for the mapping." << std::endl;
                        // THERE ARE TWO SITUATIONS WHERE WE WOULD START A NEW MAPPING...
                        // A). WE ARE MAPPING TO A COMPLETELY DIFFERENT NODE...
                        // B). THE DIRECTION OF THE MAPPING IS FLIPPED...

                        if (!mapping.mapping_continues(graph_pos)) {
                            sglib::OutputLog(sglib::LogLevels::DEBUG) << "Kmer match does not continues the current mapping." << std::endl;
                            sequence_mappings.push_back(mapping);
                            mapping.start_new_mapping(graph_pos, sk.pos, graph_kmer_index);
                        } else { // IF NODE KMER MAPS TO IS THE SAME, EXTEND THE CURRENT MAPPING...
                            sglib::OutputLog(sglib::LogLevels::DEBUG) << "Kmer match continues the current mapping." << std::endl;
                            mapping.extend(graph_pos.pos, sk.pos);
                        }
                    }
                } else {
#pragma omp critical (save_unmapped_kmer)
                    {
                        unmapped_kmers.emplace_back(sk.kmer);
                    }
                }
            }
            // TODO : Check if this last push_back is required
            if (mapping.ismatched()) sequence_mappings.push_back(mapping);

#pragma omp critical (store_seqmapping)
            {
                mappings_of_sequence[sequence.id] = sequence_mappings;
            }

#pragma omp critical (readrec)
            {
                c = fastaReader.next_record(sequence);
            }
        }
    }
    sglib::OutputLog(sglib::LogLevels::INFO) << "Mapped " << mapped_kmers_count << " Kmers from " << sequence_count << " sequences." << std::endl;
    sglib::OutputLog(sglib::LogLevels::INFO) << "Failed to map " << unmapped_kmers.size() << " unique kmers from the reference." << std::endl;
 }

void SequenceThreader::thread_mappings() {
    sglib::OutputLog(sglib::LogLevels::INFO) << "Assembling mappings into contiguous threads." << std::endl;
    for(const auto& sequence_mappings : mappings_of_sequence) {

        // For every sequence, initialize and empty sequence path, and a vector to store constructed paths.
        sglib::OutputLog(sglib::LogLevels::DEBUG) << "For sequence: " << sequence_mappings.first << std::endl;
        SequenceMappingThread mapping_thread(sg);
        std::vector<SequenceMappingThread> mapping_threads;

        for(const auto &sm : sequence_mappings.second) {

            // For every mapping hit the sequence has, try to append the node of the mapping hit to the current path.

            sglib::OutputLog(sglib::LogLevels::DEBUG) << "Considering the mapping:" << std::endl << sm << std::endl;

            bool could_append = mapping_thread.append_mapping(sm);

            if (!could_append) {
                sglib::OutputLog(sglib::LogLevels::DEBUG) << "Saving current path." << std::endl;
                mapping_threads.emplace_back(mapping_thread);
                sglib::OutputLog(sglib::LogLevels::DEBUG) << "Clearing path to start a new path." << std::endl;
                mapping_thread.clear();
                sglib::OutputLog(sglib::LogLevels::DEBUG) << "Adding node to new path." << std::endl;
                mapping_thread.append_mapping(sm);
            }
            sglib::OutputLog(sglib::LogLevels::DEBUG) << "Current path of mappings is " << mapping_thread.size() << " nodes long." << std::endl;
        }
        mapping_threads.emplace_back(mapping_thread);
        mapping_threads_of_sequence[sequence_mappings.first] = mapping_threads;
    }
}


std::vector<SequenceGraphPath> SequenceThreader::collect_paths(const sgNodeID_t seed, const sgNodeID_t target, unsigned int size_limit, unsigned int edge_limit) {

    struct visitor {
        sgNodeID_t node;
        uint dist;
        uint path_length;
        visitor(sgNodeID_t n, uint d, uint p) : node(n), dist(d), path_length(p) {}
    };

    std::vector<SequenceGraphPath> collected_paths;
    std::vector<sgNodeID_t> current_path;

    std::stack<visitor> to_visit;
    to_visit.emplace(seed, 0, 0);
    std::set<sgNodeID_t> visited;

    while (!to_visit.empty()) {
        const auto activeNode(to_visit.top());
        to_visit.pop();
        if (visited.find(activeNode.node) == visited.end() and (activeNode.path_length < edge_limit or edge_limit==0) and (activeNode.dist < size_limit or size_limit==0) )
        {
            visited.emplace(activeNode.node);

            //
            if (current_path.size() > activeNode.path_length) {
                current_path[activeNode.path_length] = activeNode.node;
                current_path.resize(activeNode.path_length + 1);
            } else {
                current_path.emplace_back(activeNode.node);
            }

            for (const auto &l : sg.get_fw_links(activeNode.node)) {
                // For each child of the current node, create a child visitor object, and enter it into the parent map.
                visitor child { l.dest, activeNode.dist + uint(sg.nodes[l.dest > 0 ? l.dest : -l.dest].sequence.length()), activeNode.path_length + 1 };

                // If I've found the target along this dfs path, I need to collect it and put it into a SequenceGraphPath object.
                if (child.node == target) {
                    current_path.emplace_back(target);
                    collected_paths.emplace_back(sg, current_path);
                } else {
                    to_visit.emplace(child);
                }
            }
        }
    }
    return collected_paths;
}

void SequenceThreader::connect_threads() {
    FastaReader<FastaRecord> fastareader({0}, query_seq_file);
    FastaRecord refsequence;
    std::vector<KmerIDX> ref_seq_kmers;
    kmerIDXFactory<FastaRecord> kmerfactory({k});
    bool c;
    c = fastareader.next_record(refsequence);
    while (c) {
        std::vector<SequenceMappingThread> threads = mapping_threads_of_sequence[refsequence.id];
        for (auto second_thread = std::next(threads.begin()); second_thread != threads.end(); second_thread++) {
            auto first_thread = std::prev(second_thread);

            auto from_mapping = first_thread->last_mapping();
            auto to_mapping = second_thread->first_mapping();

            auto first_idx = from_mapping.query_start();
            auto final_idx = to_mapping.query_end();
            if (first_idx > final_idx) {
                std::swap(first_idx, final_idx);
            }
            size_t len = size_t(final_idx) - size_t(first_idx) + 1;
            std::string refsubstr = refsequence.seq.substr((size_t) first_idx, len);

            // Try to collect some paths....
            auto graphpaths = collect_paths(from_mapping.dirnode(), to_mapping.dirnode(), 0, 10);
            if (graphpaths.size() > 0) {
                std::cout << "Found paths between two threads..." << std::endl;

                for(const auto& path : graphpaths) {
                    std::cout << path.get_fasta_header(true) << std::endl;
                }

            } else {
                std::cout << "Did not find paths between two threads..." << std::endl;
            }
        }
        c = fastareader.next_record(refsequence);
    }
}




void SequenceThreader::graph_threads_to_fasta(std::ofstream& output_file, bool use_oldnames) const {
    for (const auto& sp : mapping_threads_of_sequence) {
        for (const auto& mapping_thread : sp.second) {
            output_file << ">Sequence_" << sp.first << '[';
            output_file << mapping_thread.query_start() << ", ";
            output_file << mapping_thread.query_end() << "] ";
            mapping_thread.print_path_header(output_file, use_oldnames);
            output_file << std::endl;
            mapping_thread.print_sequence(output_file);
            output_file << std::endl;
        }
    }
}

void SequenceThreader::query_threads_to_fasta(std::ofstream& output_file, bool use_oldnames) const {
    FastaRecord sequence;
    FastaReader<FastaRecord> fastaReader({0}, query_seq_file);
    bool c;
    c = fastaReader.next_record(sequence);
    while (c) {
        const auto& sequence_threads = mapping_threads_of_sequence.at((seqID_t)sequence.id);
        for (const auto& thread : sequence_threads) {
            output_file << ">Sequence_" << sequence.id << '[';
            output_file << thread.query_start() << ", ";
            output_file << thread.query_end() << "] ";
            thread.print_path_header(output_file, use_oldnames);
            output_file << std::endl;
            uint32_t first_idx = thread.query_start();
            uint32_t final_idx = thread.query_end();
            if (first_idx > final_idx) {
                std::swap(first_idx, final_idx);
            }
            size_t len = size_t(final_idx) - size_t(first_idx) + 1;
            auto subseq = sequence.seq.substr((size_t) first_idx, len);
            output_file << subseq << std::endl;
        }
        c = fastaReader.next_record(sequence);
    }
}

void SequenceThreader::print_mappings(std::ostream& out, bool use_oldnames) const {
    for( auto it = mappings_of_sequence.begin(); it != mappings_of_sequence.end(); ++it) {
        for (const auto &sm:it->second) {
            auto sd = sm.seq_direction();
            auto nd = sm.node_direction();
            auto seqdir = sd == Forward ? "Forward" : sd == Backwards ? "Backwards" : "DIRERR";
            auto nodedir = nd == Forward ? "Forward" : nd == Backwards ? "Backwards" : "DIRERR";
            out << "Sequence: " << sm.seq_id << " from: ";
            out << sm.first_seq_pos << " : " << sm.last_seq_pos << " (" << seqdir << ")";
            out << ", maps to node: " << (use_oldnames ? sg.nodeID_to_name(sm.absnode()) : std::to_string(sm.absnode()));
            out << " from: " << sm.first_node_pos << " : " << sm.last_node_pos  << " (" << nodedir << "). ";
            out << sm.matched_unique_kmers << " / " << sm.possible_unique_matches << " unique kmers matched.";
            out << " (" << sm.n_kmers_in_node << " kmers exist in node).";
            out << std::endl;
        }
    }
}

void SequenceThreader::print_paths(std::ostream& out, bool use_oldnames) const {
    for(auto it = mapping_threads_of_sequence.begin(); it != mapping_threads_of_sequence.end(); ++it) {
        out << "Threads for query sequence: " << it->first << std::endl;
    }
}

void SequenceThreader::print_dark_nodes(std::ofstream& output_file) const {
    for(sgNodeID_t i = 1; i <= sg.nodes.size(); i++){
        if (graph_kmer_index.is_unmappable(i)) {
            output_file << sg.nodeID_to_name(i) << std::endl;
        }
    }
}

void SequenceThreader::print_full_node_diagnostics(std::ofstream &output_file) const {
    output_file << "nodeID\tseq_size\tn_kmers\tn_unique_kmers\tn_forward_links\tn_backward_links" << std::endl;
    for(sgNodeID_t i = 1; i <= sg.nodes.size(); i++){
        output_file << sg.nodeID_to_name(i) << '\t' << sg.nodes[i].sequence.size() << '\t';
        output_file << std::to_string(graph_kmer_index.total_kmers_in_node(i)) << '\t';
        output_file << std::to_string(graph_kmer_index.unique_kmers_in_node(i)) << '\t';
        output_file << std::to_string(sg.get_fw_links(i).size()) << '\t';
        output_file << std::to_string(sg.get_bw_links(i).size()) << '\t';
        output_file << std::endl;
    }
}

void SequenceThreader::print_unmapped_nodes(std::ofstream& output_file) const {
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