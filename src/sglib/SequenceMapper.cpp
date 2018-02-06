//
// Created by Ben Ward (EI) on 29/01/2018.
//

#include "SequenceMapper.h"
#include "sglib/SMR.h"
#include "sglib/factories/KMerIDXFactory.h"
#include "sglib/readers/SequenceGraphReader.h"

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
    stream << " from: " << sm.first_node_pos << " : " << sm.last_node_pos  << " (" << nodedir << ")";
    stream << ", with " << sm.unique_matches << " unique matches.";
    return stream;
}

void SequenceMapping::initiate_mapping(uint64_t sequence_id) {
    seq_id = sequence_id;
    first_seq_pos = 0;
    last_seq_pos = 0;
    node = 0;
    first_node_pos = 0;
    last_node_pos = 0;
    unique_matches = 0;
}

bool SequenceMapping::ismatched(){
    return node != 0;
}

void SequenceMapping::start_new_mapping(sgNodeID_t nodeid, int32_t nodepos, int32_t seqpos) {
    first_seq_pos = seqpos;
    last_seq_pos = seqpos;
    node = nodeid;
    first_node_pos = nodepos;
    last_node_pos = nodepos;
    unique_matches = 1;
}

void SequenceMapping::extend(int32_t nodepos, int32_t seqpos) {
    unique_matches++;
    last_node_pos = nodepos;
    last_seq_pos = seqpos;
}

sgNodeID_t SequenceMapping::absnode() const {
    return std::abs(node);
}

int32_t SequenceMapping::n_unique_matches() const {
    return unique_matches;
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

bool SequenceMapping::mapping_continues(sgNodeID_t next_node, int32_t next_position) const {
    //std::cout << "Deciding if mapping will continue..." << std::endl;
    auto same_node = absnode() == std::abs(next_node);
    //std::cout << "Same node: " << same_node << std::endl;
    auto direction_continues = direction_will_continue(next_position);
    //std::cout << "Direction continues: " << direction_continues << std::endl;
    return same_node and direction_continues;
}



void SequenceMapper::generate_unique_kmer_index(uint8_t k) {

    SMR<KmerIDX, kmerIDXFactory<FastaRecord>,
        GraphNodeReader<FastaRecord>, FastaRecord,
        GraphNodeReaderParams, KMerIDXFactoryParams> kmerIDX_SMR({1, sg}, {k}, 10*GB, 0, 1, output_prefix);

    std::cout << "Generating unique graph " << k << "mers and building index..." << std::endl;
    for (auto &kidx :kmerIDX_SMR.process_from_memory()) graph_kmer_index[kidx.kmer]={kidx.contigID,kidx.pos};
    std::vector<uint64_t> uniqKmer_statistics(kmerIDX_SMR.summaryStatistics());
    std::cout << "Number of sequences in graph: " << uniqKmer_statistics[2] << std::endl;
    std::cout << "Number of " << int(k) << "-kmers in graph " << uniqKmer_statistics[0] << std::endl;
    std::cout << "Number of " << int(k) << "-kmers in graph index " << uniqKmer_statistics[1] << std::endl;

}

void SequenceMapper::map_sequences_from_file(const uint64_t min_matches, const std::string& filename) {

    std::cout << "Mapping sequence kmers to graph..." << std::endl;
    FastaReader<FastaRecord> fastaReader({0}, filename);
    std::atomic<uint64_t> mapped_kmers_count(0), sequence_count(0);

#pragma omp parallel shared(fastaReader, mappings_in_node, mappings_of_sequence)
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
            mapping.initiate_mapping((uint64_t)sequence.id);

            // get all Kmers from sequence.
            seqkmers.clear();

            kf.setFileRecord(sequence);
            kf.next_element(seqkmers);

            // For each kmer on sequence.
            std::vector<SequenceMapping> sequence_mappings;
            for (auto &sk:seqkmers) {
                auto nk = graph_kmer_index.find(sk.kmer);
                // IF KMER EXISTS ON GRAPH
                if (graph_kmer_index.end() != nk) {
                    mapped_kmers_count++;
                    // IF THE KMER MATCH IS THE FIRST MATCH FOR THE MAPPING...
                    if (!mapping.ismatched()) {
                        mapping.start_new_mapping(nk->second.node, nk->second.pos, sk.pos);
                    } // IF THE KMER MATCH IS NOT THE FIRST MATCH FOR THE MAPPING...
                    else {
                        // THERE ARE TWO SITUATIONS WHERE WE WOULD START A NEW MAPPING...
                        // A). WE ARE MAPPING TO A COMPLETELY DIFFERENT NODE...
                        // B). THE DIRECTION OF THE MAPPING IS FLIPPED...
                        //std::cout << "Current mapping: " << mapping << std::endl;

                        if (!mapping.mapping_continues(nk->second.node, nk->second.pos)) {
                            //std::cout << "Mapping to new node, making new mapping..." << std::endl;
                            sequence_mappings.push_back(mapping);
                            mapping.start_new_mapping(nk->second.node, nk->second.pos, sk.pos);
                            //std::cout << mapping << std::endl;
                        } else { // IF NODE KMER MAPS TO IS THE SAME, EXTEND THE CURRENT MAPPING...
                            mapping.extend(nk->second.pos, sk.pos);
                        }
                    }
                }
            }
            // TODO : Check if this last push_back is required
            if (mapping.ismatched()) sequence_mappings.push_back(mapping);

            for (const auto &sm:sequence_mappings)
                if (sm.absnode() != 0 and sm.n_unique_matches() >= min_matches) {
#pragma omp critical (storeA)
                    {
                        mappings_in_node[std::abs(sm.absnode())].emplace_back(sm);
                    }
                    ++mapped_kmers_count;
                }
/*
            for (const auto &sm:sequence_mappings) {
                std::cout << sm << std::endl;
            }
*/
#pragma omp critical (storeB)
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

#pragma omp parallel for
    for (sgNodeID_t n = 1; n < mappings_in_node.size(); ++n){
        std::sort(mappings_in_node[n].begin(), mappings_in_node[n].end());
    }

    std::cout << "Mapped " << mapped_kmers_count << " Kmers from " << sequence_count << " sequences." << std::endl;
}

void SequenceMapper::mappings_paths() {
    for(auto sequence_mappings = mappings_of_sequence.begin(); sequence_mappings != mappings_of_sequence.end(); ++sequence_mappings) {
        std::cout << "Trying to make a path from mappings of sequence: " << sequence_mappings->first << std::endl;
        std::cout << "Making node list..." << std::endl;
        std::vector<sgNodeID_t> nodeList;
        for(const auto &sm:sequence_mappings->second) {
            auto dirnode = sm.node_direction() == Forward ? sm.absnode() : -sm.absnode();
            nodeList.push_back(dirnode);
        }
        std::cout << "Constructing sequence path..." << std::endl;
        SequenceGraphPath output_path(sg, nodeList);
        auto mappedseq = output_path.get_sequence();
        std::cout << "Constructed sequence!" << std::endl;
        std::cout << "Managed to map reference sequence to genome graph!" << std::endl;
        std::cout << "Writing " << mappedseq.length() << " of aligned reference to file!" << std::endl;
        std::ofstream file(output_prefix + "mapped_sequence.fasta");
        file << "> Mapped Sequence" << std::endl << mappedseq << std::endl;
        file.close();
    }

}