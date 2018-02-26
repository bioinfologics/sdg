//
// Created by Ben Ward (EI) on 29/01/2018.
//

#include <numeric>
#include <atomic>
#include <sglib/SequenceThreader.h>
#include <sglib/factories/KMerIDXFactory.h>
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
    //std::cout << "Deciding if mapping will continue..." << std::endl;
    auto same_node = absnode() == std::abs(gpos.node);
    //std::cout << "Same node: " << same_node << std::endl;
    auto direction_continues = direction_will_continue(gpos.pos);
    //std::cout << "Direction continues: " << direction_continues << std::endl;
    return same_node and direction_continues;
}

double SequenceMapping::POPUKM() const {
    return matched_unique_kmers / possible_unique_matches;
}

void SequenceThreader::map_sequences_from_file(const uint64_t min_matches, const std::string& filename) {

    sglib::OutputLog(sglib::LogLevels::INFO) << "Mapping sequence kmers to graph." << std::endl;
    FastaReader<FastaRecord> fastaReader({0}, filename);
    std::atomic<uint64_t> mapped_kmers_count(0), sequence_count(0);

#pragma omp parallel shared(fastaReader, mappings_of_sequence)
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
}

void SequenceThreader::mappings_paths() {

    sglib::OutputLog(sglib::LogLevels::INFO) << "Assembling mappings into paths." << std::endl;

    for(const auto& sequence_mappings : mappings_of_sequence) {

        // For every sequence, initialize and empty sequence path, and a vector to store constructed paths.
        sglib::OutputLog(sglib::LogLevels::DEBUG) << "For sequence: " << sequence_mappings.first << std::endl;
        SequenceGraphPath sgpath(sg);
        std::vector<std::vector<SequenceMapping>> seq_mapping_paths(0);
        std::vector<SequenceMapping> mapping_path(0);

        for(const auto &sm:sequence_mappings.second) {

            // For every mapping hit the sequence has, try to append the node of the mapping hit to the current path.

            sglib::OutputLog(sglib::LogLevels::DEBUG) << "Considering the mapping:" << std::endl << sm << std::endl;

            auto dn = sm.dirnode();

            sglib::OutputLog(sglib::LogLevels::DEBUG) << "Trying to add to path as: " << dn << std::endl;

            bool could_append = sgpath.append_to_path(dn);

            if (could_append) {
                sglib::OutputLog(sglib::LogLevels::DEBUG) << "Was able to append " << dn << " to the path." << std::endl;
                mapping_path.emplace_back(sm);
            } else {
                sglib::OutputLog(sglib::LogLevels::DEBUG) << "Was not able to append " << dn << " to the path." << std::endl;
                sglib::OutputLog(sglib::LogLevels::DEBUG) << "Saving current path." << std::endl;
                seq_mapping_paths.emplace_back(mapping_path);
                sglib::OutputLog(sglib::LogLevels::DEBUG) << "Clearing path to start a new path." << std::endl;
                mapping_path.clear();
                sgpath.nodes.clear();
                sglib::OutputLog(sglib::LogLevels::DEBUG) << "Adding node to new path." << std::endl;
                sgpath.append_to_path(dn);
                mapping_path.emplace_back(sm);
            }
            sglib::OutputLog(sglib::LogLevels::DEBUG) << "Current path of mappings is " << mapping_path.size() << " nodes long." << std::endl;
        }

        seq_mapping_paths.emplace_back(mapping_path);
        paths_of_mappings_of_sequence[sequence_mappings.first] = seq_mapping_paths;

    }

}

void SequenceThreader::graph_paths_to_fasta(std::ofstream& output_file) const {
    for (const auto& sp : paths_of_mappings_of_sequence) {
        for (const auto& mapping_path : sp.second) {
            std::vector<sgNodeID_t> nodes;
            SequenceGraphPath p(sg);
            output_file << ">Sequence_" << sp.first << "_path";
            for (const auto& sm : mapping_path) {
                auto dn = sm.dirnode();
                p.nodes.emplace_back(dn);
                char dir = dn > 0 ? '+' : '-';
                output_file << ' ' << dir << sg.nodeID_to_name(std::abs(dn));
            }
            auto sequence = p.get_sequence();
            output_file << std::endl << sequence << std::endl;
        }
    }
}

void SequenceThreader::print_mapping_path_name(const std::vector<SequenceMapping>& path, std::ofstream& output_file) const {
    output_file << "[ ";
    for (const auto& mapping : path) {
        output_file << mapping.dirnode() << ", ";
    }
    output_file << " ]";
}

void SequenceThreader::query_paths_to_fasta(std::ofstream& output_file) const {
    FastaRecord sequence;
    FastaReader<FastaRecord> fastaReader({0}, query_seq_file);
    bool c;
    c = fastaReader.next_record(sequence);
    while (c) {
        const std::vector<std::vector<SequenceMapping>>& mapping_paths = paths_of_mappings_of_sequence.at((seqID_t)sequence.id);
        for (const auto& mapping_path : mapping_paths) {
            auto first_idx = mapping_path.front().first_seq_pos;
            auto final_idx = mapping_path.back().last_seq_pos;
            output_file << sequence.name << '[' << first_idx << " -> " << final_idx << "] ";
            print_mapping_path_name(mapping_path, output_file);
            if (first_idx > final_idx) std::swap(first_idx, final_idx);
            size_t len = size_t(final_idx) - size_t(first_idx) + 1;
            auto subseq = sequence.seq.substr((size_t) first_idx, len);
            output_file << std::endl << subseq << std::endl;
        }
        c = fastaReader.next_record(sequence);
    }
}

std::set<SequenceGraphPath> SequenceThreader::all_unique_paths() const {
    std::set<SequenceGraphPath> paths;
    for(const auto& sp : paths_of_mappings_of_sequence) {
        for(const auto& mapping_path : sp.second) {
            std::vector<sgNodeID_t> nodes;
            SequenceGraphPath sgpath(sg);
            for(const auto& sm : mapping_path) {
                auto dn = sm.dirnode();
                sgpath.nodes.emplace_back(dn);
            }
            paths.emplace(sgpath);
        }
    }
    return paths;
}

void SequenceThreader::print_unique_paths_sizes(std::ofstream& output_file) const {
    std::set<SequenceGraphPath> unique_paths = all_unique_paths();
    for(auto& path : unique_paths) {
        std::string seq = path.get_sequence();
        auto seqsize = seq.length();
        auto pathname = path.get_fasta_header().erase(0, 1);
        output_file << pathname << '\t' << seqsize << std::endl;
    }
}

void SequenceThreader::print_mappings() const {
    for( auto it = mappings_of_sequence.begin(); it != mappings_of_sequence.end(); ++it) {
        sglib::OutputLog(sglib::LogLevels::DEBUG) << "Mappings for query sequence: " << it->first << ":" << std::endl;
        if(sglib::OutputLogLevel == sglib::LogLevels::DEBUG)
            for (const auto &sm:it->second) {
                std::cout << sm << std::endl;
            }
    }
}

void SequenceThreader::print_paths() const {
    for( auto it = paths_of_mappings_of_sequence.begin(); it != paths_of_mappings_of_sequence.end(); ++it) {
        sglib::OutputLog(sglib::LogLevels::DEBUG) << "Paths for query sequence: " << it->first << std::endl;
        if(sglib::OutputLogLevel == sglib::LogLevels::DEBUG)
            for (const auto &path:it->second) {
                std::cout << "Path of " << path.size() << " mappings: ";
                for (const auto &sm:path) {
                    sgNodeID_t dirnode = sm.node_direction() == Forward ? sm.absnode() : -sm.absnode();
                    std::cout << dirnode << ", ";
                }
                std::cout << std::endl;
            }
    }
}