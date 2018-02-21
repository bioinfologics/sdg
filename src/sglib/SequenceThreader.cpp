//
// Created by Ben Ward (EI) on 29/01/2018.
//

#include "SequenceThreader.h"
#include "factories/KMerIDXFactory.h"
#include <numeric>
#include <atomic>

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

    std::cout << "Mapping sequence kmers to graph..." << std::endl;
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
                    std::cout << '(' << graph_pos.node << ", " << graph_pos.pos << "), ";
                    mapped_kmers_count++;
                    // IF THE KMER MATCH IS THE FIRST MATCH FOR THE MAPPING...
                    if (!mapping.ismatched()) {
                        mapping.start_new_mapping(graph_pos, sk.pos, graph_kmer_index);
                    } // IF THE KMER MATCH IS NOT THE FIRST MATCH FOR THE MAPPING...
                    else {
                        // THERE ARE TWO SITUATIONS WHERE WE WOULD START A NEW MAPPING...
                        // A). WE ARE MAPPING TO A COMPLETELY DIFFERENT NODE...
                        // B). THE DIRECTION OF THE MAPPING IS FLIPPED...

                        if (!mapping.mapping_continues(graph_pos)) {
                            sequence_mappings.push_back(mapping);
                            mapping.start_new_mapping(graph_pos, sk.pos, graph_kmer_index);
                        } else { // IF NODE KMER MAPS TO IS THE SAME, EXTEND THE CURRENT MAPPING...
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

            auto tc = ++sequence_count;
            if (tc % 100000 == 0) std::cout << mapped_kmers_count << " / " << tc << std::endl;

            #pragma omp critical (readrec)
            {
                c = fastaReader.next_record(sequence);
            }
        }
    }

    std::cout << "Mapped " << mapped_kmers_count << " Kmers from " << sequence_count << " sequences." << std::endl;
}

void SequenceThreader::mappings_paths() {
    for(const auto& sequence_mappings : mappings_of_sequence) {

        // For every sequence, initialize and empty sequence path, and a vector to store constructed paths.
        std::cout << "Constructing mapping paths for sequence: " << sequence_mappings.first << std::endl;
        SequenceGraphPath sgpath(sg);
        std::vector<std::vector<SequenceMapping>> seq_mapping_paths(0);
        std::vector<SequenceMapping> mapping_path(0);

        for(const auto &sm:sequence_mappings.second) {

            // For every mapping hit the sequence has, try to append the node of the mapping hit to the current path.

            //std::cout << "CONSIDERING THE FOLLOWING MAPPING:" << std::endl << sm << std::endl;

            //sgNodeID_t dirnode = sm.node_direction() == Forward ? sm.absnode() : -sm.absnode();
            auto dn = sm.dirnode();

            //std::cout << "Trying to add to path as: " << dn << std::endl;

            bool could_append = sgpath.append_to_path(dn);

            if (could_append) {
                //std::cout << "Was able to append " << dn << " to the path." << std::endl;
                //std::cout << "Adding mapping to the path of mappings." << std::endl;
                mapping_path.emplace_back(sm);
                //std::cout << "Current path of mappings is " << mapping_path.size() << " nodes long." << std::endl;
            } else {
                //std::cout << "Was not able to append " << dn << " to the path." << std::endl;
                //std::cout << "Saving current path." << std::endl;
                seq_mapping_paths.emplace_back(mapping_path);
                //std::cout << "Clearing path to start a new path." << std::endl;
                mapping_path.clear();
                sgpath.nodes.clear();
                //std::cout << "ADDING NODE TO NEW PATH..." << std::endl;
                sgpath.append_to_path(dn);
                mapping_path.emplace_back(sm);
                //std::cout << "Current path of mappings is " << mapping_path.size() << " nodes long." << std::endl;
            }
        }

        seq_mapping_paths.emplace_back(mapping_path);
        paths_of_mappings_of_sequence[sequence_mappings.first] = seq_mapping_paths;

    }

}

void SequenceThreader::paths_to_fasta(std::ofstream& output_file) const {
    for(const auto& sp : paths_of_mappings_of_sequence) {
        for(const auto& mapping_path : sp.second) {
            std::vector<sgNodeID_t> nodes;
            SequenceGraphPath p(sg);
            output_file << ">Sequence_" << sp.first << "_path";
            for(const auto& sm : mapping_path) {
                auto dn = sm.dirnode();
                p.nodes.emplace_back(dn);
                output_file << '_' << dn;
            }
            auto sequence = p.get_sequence();
            output_file << std::endl << sequence << std::endl;
        }
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