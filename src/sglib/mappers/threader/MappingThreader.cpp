//
// Created by Ben Ward (EI) on 21/05/2018.
//

#include <sglib/mappers/threader/MappingThreader.h>
#include <sglib/processors/PathExplorer.h>
#include <sglib/utilities/nstats.h>

// Mapping Threader...

// Public methods...

void MappingThreader::thread_mappings() {
    sglib::OutputLog() << "Assembling mappings into contiguous threads" << std::endl;

    FastaReader<FastaRecord> fastareader({0}, getQueryFile());
    FastaRecord refsequence;
    std::vector<KmerIDX> ref_seq_kmers;
    kmerIDXFactory<FastaRecord> kmerfactory({ getK() });
    PathExplorer pe { getGraph() };
    MappingThread mappingThread(getGraph());
    std::vector<MappingThread> mappingThreads;

    const std::unordered_map<seqID_t, std::vector<NodeMapping>> &allMappings { mapper.getMappings() };

    bool c;
    c = fastareader.next_record(refsequence);

    while (c) {
        sglib::OutputLog() << "Threading mappings for " << refsequence.name << std::endl;
        const auto &mappings = allMappings.at(refsequence.id);

        for(const auto &sm : mappings) {
            // For every mapping hit the sequence has, try to append the node of the mapping hit to the current path.
            std::cout << "Considering the mapping:" << std::endl << sm << std::endl;

            std::cout << "Can this mapping be connected to the previous one trivially?" << std::endl;

            // There are a few cases that could happen here. First, either the mapping could be appended, trivially.

            bool could_append_trivially = append_mapping_trivially(mappingThread, sm);

            std::cout << (could_append_trivially ? " Yes" : " No") << std::endl;

            if(!could_append_trivially) {
                std::cout << "Can this mapping be connected through path finding?" << std::endl;

                // Work out the region of the sub sequence to extract from the reference sequence...
                auto first_idx = mappingThread.lastMapping().query_start();
                auto final_idx = sm.query_end();
                if (first_idx > final_idx) {
                    std::swap(first_idx, final_idx);
                }
                // Extract the sub-sequence from the reference sequence...
                size_t len { size_t(final_idx) - size_t(first_idx) + 1 };
                std::string refsubstr { refsequence.seq.substr((size_t) first_idx, len) };

                SequenceGraphPath bestpath { getGraph() };
                bool found_path = pe.find_best_path(bestpath, mappingThread.lastMapping().dirnode(), sm.dirnode(), refsubstr);

                if (found_path) {
                    std::cout << "Found a suitable path connecting the two mappings" << std::endl;
                    std::cout << "Adding the nodes of this suitable path to the sequence mapping thread" << std::endl;
                    std::cout << "Current bestpath: ";
                    for (const auto &node : bestpath.nodes) {
                        std::cout << node << ", ";
                    }
                    std::cout << std::endl;
                    std::cout << "Adding nodes of bestpath to the current thread:" << std::endl;
                    append_mapping_with_path(mappingThread, sm, bestpath);
                } else {
                    // What happens if you couldn't append trivially, AND could not append by searching a path.
                    std::cout << "Could not trivially link mappings, or get there through path search" << std::endl;
                    sglib::OutputLog(sglib::LogLevels::DEBUG) << "Saving current path." << std::endl;
                    mappingThreads.emplace_back(mappingThread);
                    sglib::OutputLog(sglib::LogLevels::DEBUG) << "Clearing path to start a new path." << std::endl;
                    mappingThread.clear();
                    sglib::OutputLog(sglib::LogLevels::DEBUG) << "Adding node to new path." << std::endl;
                    append_mapping_trivially(mappingThread, sm);
                }
            }
            std::cout << "Current path of mappings is " << mappingThread.size() << " nodes long." << std::endl;
        }
        mapping_threads_of_sequence.emplace(refsequence.id, mappingThreads);
        c = fastareader.next_record(refsequence);
    }
}

void MappingThreader::graph_threads_to_fasta(std::ofstream& output_file, bool use_oldnames) const {
    for (const auto& sp : mapping_threads_of_sequence) {
        for (const auto& mapping_thread : sp.second) {
            output_file << ">Sequence_" << sp.first << '[';
            output_file << mapping_thread.queryStart() << ", ";
            output_file << mapping_thread.queryEnd() << "] ";
            mapping_thread.print_path_header(output_file, use_oldnames);
            output_file << std::endl;
            mapping_thread.print_sequence(output_file);
            output_file << std::endl;
        }
    }
}

void MappingThreader::query_threads_to_fasta(std::ofstream& output_file, bool use_oldnames) const {
    FastaRecord sequence;
    FastaReader<FastaRecord> fastaReader({0}, getQueryFile());
    bool c;
    c = fastaReader.next_record(sequence);
    while (c) {
        const auto& sequence_threads = mapping_threads_of_sequence.at((seqID_t)sequence.id);
        for (const auto& thread : sequence_threads) {
            output_file << ">Sequence_" << sequence.id << '[';
            output_file << thread.queryStart() << ", ";
            output_file << thread.queryEnd() << "] ";
            thread.print_path_header(output_file, use_oldnames);
            output_file << std::endl;
            uint32_t first_idx = thread.queryStart();
            uint32_t final_idx = thread.queryEnd();
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

// Private methods...

bool MappingThreader::append_mapping_trivially(MappingThread &thread, NodeMapping mapping) {
    auto dn = mapping.dirnode();
    sglib::OutputLog(sglib::LogLevels::DEBUG) << "Trying to add to path trivially as: " << dn << std::endl;
    auto pathSuccess = thread.nodePath().append_to_path(dn);
    if (pathSuccess) {
        sglib::OutputLog(sglib::LogLevels::DEBUG) << "Was able to trivially append " << dn << " to the path." << std::endl;
        thread.orderedMappings().emplace_back(mapping);
    } else {
        sglib::OutputLog(sglib::LogLevels::DEBUG) << "Was not able to trivially append " << dn << " to the path." << std::endl;
    }
    return pathSuccess;
}

bool MappingThreader::append_mapping_with_path(MappingThread &thread, NodeMapping mapping, SequenceGraphPath &path) {
    for (auto node = ++path.nodes.begin(); node != path.nodes.end(); ++node) {
        std::cout << "Adding node: " << *node << std::endl;
        auto status = thread.node_path.append_to_path(*node);
        std::cout << status << std::endl;
        if (!status) {
            throw std::logic_error("Appending a node during best path stuff");
        }
    }
    thread.orderedMappings().emplace_back(mapping);
}

void MappingThreader::calculate_reference_inclusion() {
    const auto& sequenceNames { mapper.querySeqNames() };
    const auto& sequenceSizes { mapper.querySeqSizes() };

    std::cout << std::endl << "///// Results Per Query Sequence /////" << std::endl;
    for (const auto& mapping_threads : mapping_threads_of_sequence) {
        std::cout << mapping_threads.first << std::endl;
        std::cout << "// " << sequenceNames.at(mapping_threads.first) << ":" << std::endl;
        std::vector<uint32_t> sizes;
        unsigned long long total_inclusion = 0;
        uint32_t current_size;
        for (const auto& mapping_thread : mapping_threads.second) {
            current_size = mapping_thread.querySize();
            total_inclusion += current_size;
            sizes.emplace_back(current_size);
        }

        auto N50 = NX(sizes, 50);
        auto N90 = NX(sizes, 90);

        auto seqSize { sequenceSizes.at(mapping_threads.first) };

        std::cout << total_inclusion << "bp / " << seqSize << "bp";
        std::cout << " ( " << ((long double) total_inclusion / seqSize) * 100;
        std::cout << "% )" << " could be threaded through the graph:" << std::endl;
        std::cout << "N50: " << N50 << "bp, N90: " << N90 << "bp." << std::endl << std::endl;
    }
}

// MappingThread...
/*
bool MappingThread::append_mapping_trivially(NodeMapping mapping) {
    auto dn = mapping.dirnode();
    sglib::OutputLog(sglib::LogLevels::DEBUG) << "Trying to add to path trivially as: " << dn << std::endl;
    auto path_success = node_path.append_to_path(dn);
    if (path_success) {
        sglib::OutputLog(sglib::LogLevels::DEBUG) << "Was able to trivially append " << dn << " to the path." << std::endl;
        ordered_mappings.emplace_back(mapping);
    } else {
        sglib::OutputLog(sglib::LogLevels::DEBUG) << "Was not able to trivially append " << dn << " to the path." << std::endl;
    }
    return path_success;
}

bool MappingThread::append_mapping_with_path(NodeMapping mapping, SequenceGraphPath &path) {
    for (auto node = ++path.nodes.begin(); node != path.nodes.end(); ++node) {
        std::cout << "Adding node: " << *node << std::endl;
        auto status = append_node(*node);
        std::cout << status << std::endl;
        if (!status) {
            throw std::logic_error("Appending a node during best path stuff");
        }
    }
    ordered_mappings.emplace_back(mapping);
}
*/

uint32_t MappingThread::querySize() const {
    const auto first = queryStart();
    const auto second = queryEnd();
    const auto sub = second > first ? second - first : first - second;
    return sub + 1;
}

void MappingThread::print_path_header(std::ostream& output_file, bool use_oldnames) const {
    auto fasta_header = node_path.get_fasta_header(use_oldnames);
    fasta_header.erase(fasta_header.begin());
    output_file << fasta_header;
}

void MappingThread::print_sequence(std::ofstream& output_file) const {
    auto s = node_path.get_sequence();
    output_file << s;
}