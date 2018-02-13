//
// Created by Luis Yanes (EI) on 12/02/2018.
//

#ifndef BSG_LONGREADMAPPER_H
#define BSG_LONGREADMAPPER_H


#include <iostream>
#include <sglib/factories/KMerIDXFactory.h>
#include <sglib/readers/FileReader.h>
#include <sglib/readers/SequenceGraphReader.h>
#include <sglib/SMR.h>
#include <sglib/PairedReadMapper.h>

class LongReadMapper {
    SequenceGraph & sg;
    uint8_t k=31;
    std::unordered_map<uint64_t, graphPosition> kmer_to_graphposition;
    std::vector<std::vector<sgNodeID_t>> nodes_in_read;    // Nodes in read
    std::vector<std::vector<size_t >> reads_in_node;      // Reads in node
public:
    LongReadMapper(uint k, SequenceGraph sg) : sg(sg), k(k) {
        update_graph_index();
    }

    void update_graph_index() {
        const int max_coverage = 1;
        const std::string output_prefix("./");
        SMR<KmerIDX,
                kmerIDXFactory<FastaRecord>,
                GraphNodeReader<FastaRecord>,
                FastaRecord, GraphNodeReaderParams, KMerIDXFactoryParams> kmerIDX_SMR({1, sg}, {k}, {4*GB, 0, max_coverage,
                                                                                                     output_prefix});

        // Get the unique_kmers from the graph into a map
        std::cout << "Indexing graph... " << std::endl;
        auto kmers(kmerIDX_SMR.process_from_memory());

        for (auto &kidx : kmers) kmer_to_graphposition[kidx.kmer]={kidx.contigID,kidx.pos};

        std::vector<uint64_t> uniqKmer_statistics(kmerIDX_SMR.summaryStatistics());
        std::cout << "Number of sequences in graph: " << uniqKmer_statistics[2] << std::endl;
        std::cout << "Number of " << int(k) << "-kmers in graph " << uniqKmer_statistics[0] << std::endl;
        std::cout << "Number of " << int(k) << "-kmers in graph index " << uniqKmer_statistics[1] << std::endl;
    }

    uint64_t map_reads(const uint16_t min_matches, std::string filename) {
        if (min_matches < 2 ) {
            throw std::invalid_argument("min_matches < 2, long reads can only be considered mapped if > 2 kmers are shared between read and sequence");
        }
        /*
         * LongRead mapping in parallel
         */
        FastqReader<FastqRecord> fastqReader({0},filename);
        std::atomic<uint64_t> mapped_count(1);
        std::atomic<uint64_t> total_count(0);
        unsigned int kmers_found=0;
#pragma omp parallel shared(fastqReader) reduction(+:kmers_found)
        {
            FastqRecord read;
            std::vector<KmerIDX> readkmers;
            kmerIDXFactory<FastqRecord> kf({this->k});
            ReadMapping mapping;
            bool c ;
#pragma omp critical(lr_fastq_reader)
            {
                c = fastqReader.next_record(read);
            }
            while (c) {
                mapping.read_id = read.id;
                //this enables partial read re-mapping by setting read_to_node to 0
                if (nodes_in_read.size()<=mapping.read_id or nodes_in_read[mapping.read_id].empty()) {
                    mapping.node = 0;
                    mapping.unique_matches = 0;
                    mapping.first_pos = 0;
                    mapping.last_pos = 0;
                    mapping.rev = false;
                    mapping.unique_matches = 0;
                    //get all kmers from read
                    readkmers.clear();

                    kf.setFileRecord(read);
                    kf.next_element(readkmers);

                    // For each kmer on read
                    std::vector<ReadMapping> read_mappings;
                    for (auto &rk:readkmers) {
                        const auto nk = kmer_to_graphposition.find(rk.kmer);
                        // If kmer exists on graph
                        if (nk != kmer_to_graphposition.cend()) {
                            kmers_found++;
                            // If first match
                            if (mapping.node == 0) {
                                mapping.node = nk->second.node;
                                if ((nk->second.node > 0 and rk.contigID > 0) or
                                    (nk->second.node < 0 and rk.contigID < 0))
                                    mapping.rev = false;
                                else mapping.rev = true;
                                mapping.first_pos = nk->second.pos;
                                mapping.last_pos = nk->second.pos;
                                mapping.unique_matches = 1;
                            }// If not first match
                            else {
                                // If node is different, end match... Start new one
                                if (mapping.node != nk->second.node) {
                                    read_mappings.push_back(mapping);
                                    mapping.node = nk->second.node;
                                    mapping.first_pos = nk->second.pos;
                                    mapping.last_pos = nk->second.pos;
                                    mapping.unique_matches=1;
                                } else {
                                    mapping.unique_matches++;
                                    mapping.last_pos = nk->second.pos;
                                }
                            }
                        }
                    }
                    // TODO : Check if this last push_back is required
                    if (mapping.node != 0) read_mappings.push_back(mapping);

                    for (auto &rm:read_mappings)
                        if (rm.node != 0 and rm.unique_matches >= min_matches) {
#pragma omp critical(lr_reads_in_node)
                            {
                                nodes_in_read[read.id].push_back(rm.node);
                            }
                            mapped_count++;
                        }
                }
                auto tc = ++total_count;
                if (tc % 100000 == 0) std::cout << mapped_count << " / " << tc << std::endl;
#pragma omp critical(lr_fastq_reader)
                {
                    c = fastqReader.next_record(read);
                }
            }

        }
        std::cout << k << "-mers found " << kmers_found <<std::endl;
        std::cout<<"Reads mapped: "<<mapped_count<<" / "<<total_count<<std::endl;

        for (auto read = nodes_in_read.cbegin(); read != nodes_in_read.cend(); ++read){
            for (const auto &node:*read) {
                reads_in_node[node].push_back(read - nodes_in_read.cbegin());
            }
        }

        return mapped_count;
    }
};


#endif //BSG_LONGREADMAPPER_H
