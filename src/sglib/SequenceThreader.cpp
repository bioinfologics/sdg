//
// Created by Ben Ward (EI) on 29/01/2018.
//

#include <numeric>
#include <atomic>
#include <algorithm>
#include <stack>
#include <vector>
#include <sglib/SequenceThreader.h>
#include <sglib/readers/FileReader.h>
#include <sglib/logger/OutputLog.h>
#include <sglib/aligners/algorithms/SmithWaterman.h>
#include <sglib/processors/PathExplorer.h>
#include <sglib/utilities/nstats.h>

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

void SequenceMapping::start_new_mapping(const graphPosition& gpos, uint32_t seqpos, const UniqueKmerIndex& index) {
    first_seq_pos = seqpos;
    last_seq_pos = seqpos;
    node = gpos.node;
    first_node_pos = gpos.pos;
    last_node_pos = gpos.pos;
    matched_unique_kmers = 1;
    possible_unique_matches = index.unique_kmers_in_node(gpos.node);
    n_kmers_in_node = index.total_kmers_in_node(gpos.node);
}

void SequenceMapping::extend(int32_t nodepos, uint32_t seqpos) {
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

#pragma omp parallel shared(fastaReader, mappings_of_sequence, unmapped_kmers, query_seq_sizes, query_seq_names)
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
                query_seq_sizes[sequence.id] = sequence.seq.size();
                query_seq_names[sequence.id] = sequence.name;
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

void SequenceThreader::filter_mappings(double score) {
    sglib::OutputLog(sglib::LogLevels::INFO) << "Filtering mappings using a threshold of " << score << "%." << std::endl;
    for (auto& sequence_mappings : mappings_of_sequence) {
        auto start { sequence_mappings.second.begin() };
        auto end { sequence_mappings.second.end() };
        auto newend {std::remove_if(start, end, [score](SequenceMapping sm) -> bool { return sm.match_score() < score; })};
        sequence_mappings.second.erase(newend, sequence_mappings.second.end());
    }
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

void SequenceThreader::bridge_threads() {

    sglib::OutputLog(sglib::LogLevels::INFO) << "Connecting fragmented mapping threads into longer paths." << std::endl;

    FastaReader<FastaRecord> fastareader({0}, query_seq_file);
    FastaRecord refsequence;
    std::vector<KmerIDX> ref_seq_kmers;
    kmerIDXFactory<FastaRecord> kmerfactory({k});
    bool c;
    c = fastareader.next_record(refsequence);

    while (c) {

        std::vector<SequenceMappingThread> threads = mapping_threads_of_sequence[refsequence.id];
        std::vector<BridgedMappingThreads> bridged_threads;

        BridgedMappingThreads growing_thread { sg, threads[0] };

        int iteration { 0 };

        auto mapping_thread { threads.begin() };
        ++mapping_thread;
        for (; mapping_thread != threads.end(); ++mapping_thread) {
            iteration++;
            auto from_mapping { growing_thread.last_mapping() };
            auto to_mapping { mapping_thread -> first_mapping() };

            // Work out the region of the sub sequence to extract from the reference sequence...
            auto first_idx = from_mapping.query_start();
            auto final_idx = to_mapping.query_end();
            if (first_idx > final_idx) {
                std::swap(first_idx, final_idx);
            }

            // Extract the subsequence from the reference sequence...
            size_t len { size_t(final_idx) - size_t(first_idx) + 1 };
            std::string refsubstr { refsequence.seq.substr((size_t) first_idx, len) };

            PathExplorer pe {sg};
            SequenceGraphPath bestpath {sg};
            int findstatus = pe.find_best_path(bestpath, from_mapping.dirnode(), to_mapping.dirnode(), refsubstr);

            // If findstatus != 0, then there was a failiure to find a path between the two nodes.
            // If bridgestatus != 0, then there is a failiure to connect paths up - likely a bug.

            int bridgestatus = 1;

            if(findstatus == 0) {
                bridgestatus = growing_thread.bridge_to_thread(bestpath, *mapping_thread);
            }

            if (bridgestatus > 0) {
                bridged_threads.emplace_back(growing_thread);
                growing_thread = { sg, *mapping_thread };
            }
        }

        bridged_mapping_threads_of_sequence[refsequence.id] = bridged_threads;

        c = fastareader.next_record(refsequence);
    }
    sglib::OutputLog(sglib::LogLevels::INFO) << "Finished connecting fragmented mapping threads." << std::endl;
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

void SequenceThreader::bridged_graph_threads_to_fasta(std::ofstream& output_file, bool use_oldnames) const {
    for (const auto& bsp : bridged_mapping_threads_of_sequence) {
        for (const auto& bridged_mapping_thread : bsp.second) {
            output_file << ">Sequence_" << bsp.first << '[';
            output_file << bridged_mapping_thread.query_start() << ", ";
            output_file << bridged_mapping_thread.query_end() << "] ";
            const auto path { bridged_mapping_thread.get_complete_path() };
            auto hdr { path.get_fasta_header(use_oldnames) };
            hdr.erase(hdr.begin());
            output_file << hdr;
            output_file << std::endl;
            output_file << path.get_sequence();
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
    for(auto it = mappings_of_sequence.begin(); it != mappings_of_sequence.end(); ++it) {
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

void SequenceThreader::calculate_reference_inclusion() {
    std::cout << std::endl << "///// Results Per Query Sequence /////" << std::endl;
    for (const auto& bridged_thread_store : bridged_mapping_threads_of_sequence) {

        std::cout << "// " << query_seq_names[bridged_thread_store.first] << ":" << std::endl;

        std::vector<uint32_t> sizes;
        unsigned long long total_inclusion = 0;
        uint32_t current_size;

        for (const auto& bridged_thread : bridged_thread_store.second) {
            current_size = bridged_thread.query_size();
            total_inclusion += current_size;
            sizes.emplace_back(current_size);
        }

        auto N50 = NX(sizes, 50);
        auto N90 = NX(sizes, 90);

        std::cout << total_inclusion << "bp / " << query_seq_sizes[bridged_thread_store.first] << "bp";
        std::cout << " ( " << ((long double) total_inclusion / query_seq_sizes[bridged_thread_store.first]) * 100;
        std::cout << "% )" << " could be threaded through the graph:" << std::endl;
        std::cout << "N50: " << N50 << "bp, N90: " << N90 << "bp." << std::endl << std::endl;
    }
}

// SequenceMappingThread...

bool SequenceMappingThread::append_mapping(SequenceMapping mapping) {
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

void SequenceMappingThread::clear() {
    ordered_mappings.clear();
    node_path.clear();
}

size_t SequenceMappingThread::size() const {
    return ordered_mappings.size();
}

uint32_t SequenceMappingThread::query_start() const {
    return ordered_mappings.front().query_start();
}

uint32_t SequenceMappingThread::query_end() const {
    return ordered_mappings.back().query_end();
}

void SequenceMappingThread::print_path_header(std::ostream& output_file, bool use_oldnames) const {
    auto fasta_header = node_path.get_fasta_header(use_oldnames);
    fasta_header.erase(fasta_header.begin());
    output_file << fasta_header;
}

void SequenceMappingThread::print_sequence(std::ofstream& output_file) const {
    auto s = node_path.get_sequence();
    output_file << s;
}

SequenceMapping SequenceMappingThread::first_mapping() const {
    return ordered_mappings.front();
}

SequenceMapping SequenceMappingThread::last_mapping() const {
    return ordered_mappings.back();
}


// Bridged Mapping Threads...

BridgedMappingThreads& BridgedMappingThreads::operator=(const BridgedMappingThreads other) {
    if (&other == this) {
        return *this;
    }
    mapping_threads = other.mapping_threads;
    bridging_paths = other.bridging_paths;
    sg = other.sg;
    return *this;
}

SequenceMapping BridgedMappingThreads::last_mapping()
{
    return mapping_threads.back().last_mapping();
}

uint32_t BridgedMappingThreads::query_start() const {
    return mapping_threads.front().query_start();
}

uint32_t BridgedMappingThreads::query_end() const {
    return mapping_threads.back().query_end();
}

int BridgedMappingThreads::bridge_to_thread(const SequenceGraphPath& sgp, const SequenceMappingThread& smt)
{
    if (mapping_threads.empty()) {
        std::cout << "Checking mapping bridge:" << std::endl;
        std::cout << sgp.get_fasta_header(true) << std::endl;
        std::cout << "ERR: No mapping thread to bridge from...";
        return 2;
    }

    const auto last_mthread_node = last_mapping().dirnode();
    if (last_mthread_node != sgp.nodes.front()) {
        std::cout << "Checking mapping bridge:" << std::endl;
        std::cout << sgp.get_fasta_header(true) << std::endl;
        std::cout << "ERR: First node of bridge mismatches last node of last mapping thread..." << std::endl;
        std::cout << last_mthread_node << " != " << sgp.nodes.front() << std::endl;
        return 3;
    }

    if (sgp.nodes.back() != smt.first_mapping().dirnode()) {
        std::cout << "Checking mapping bridge:" << std::endl;
        std::cout << sgp.get_fasta_header(true) << std::endl;
        std::cout << "ERR: First node of new mapping thread mismatches last node of bridge..." << std::endl;
        std::cout << sgp.nodes.back() << " != " << smt.first_mapping().dirnode() << std::endl;
        return 4;
    }

    bridging_paths.emplace_back(sgp);
    mapping_threads.emplace_back(smt);
    return 0;
}

SequenceGraphPath BridgedMappingThreads::get_complete_path() const {
    SequenceGraphPath finalpath { sg };

    auto current_thread { mapping_threads.cbegin() };
    auto current_thread_path { current_thread -> get_graph_path() };

    finalpath.nodes.insert(finalpath.nodes.end(), current_thread_path.nodes.cbegin(), current_thread_path.nodes.cend());

    for (const auto& bridge : bridging_paths) {
        finalpath.nodes.insert(finalpath.nodes.end(), std::next(bridge.nodes.cbegin()), std::prev(bridge.nodes.cend()));
        std::advance(current_thread, 1);
        current_thread_path = current_thread -> get_graph_path();
        finalpath.nodes.insert(finalpath.nodes.end(), current_thread_path.nodes.cbegin(), current_thread_path.nodes.cend());
    }
    return finalpath;
}