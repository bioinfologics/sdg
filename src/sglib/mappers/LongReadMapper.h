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
#include <sglib/factories/StrandedMinimiserSketchFactory.h>

class LongReadMapper {

    struct Hit {
        sgNodeID_t node;
        bool same_strand;
        int32_t intercept;
        int32_t pos;
        Hit(sgNodeID_t node, bool strand, int32_t intercept, int32_t pos) : node(node), same_strand(strand), intercept(intercept), pos(pos) {}

        bool operator<(const Hit &o) const {
            return std::tie(node,same_strand,intercept,pos) < std::tie(o.node,o.same_strand,o.intercept,o.pos);
        }
        friend std::ostream& operator<<(std::ostream& os, const Hit& hit) {
            os << "("<< hit.node << ", " << hit.intercept << ", " << hit.pos << ")";
            return os;
        }

    };

    struct graphStrandPos{
        sgNodeID_t node;
        int32_t pos;

        graphStrandPos(sgNodeID_t node, int32_t pos) : node(node), pos(pos) {}
    };

    SequenceGraph & sg;
    uint8_t k=31;
    uint8_t w=5;
    std::unordered_map<uint64_t, std::vector<graphStrandPos>> kmer_to_graphposition;
public:
    std::vector<std::vector<sgNodeID_t>> nodes_in_read;    // Nodes in read
    std::vector<std::vector<size_t >> reads_in_node;      // Reads in node

    LongReadMapper(uint8_t k, uint8_t w, SequenceGraph &sg) : sg(sg), k(k), w(w) {
        build_minimiser_index();
    }

    void build_minimiser_index() {
        StrandedMinimiserSketchFactory kf(k, w);
        std::set<MinPosIDX> sketch;
        GraphNodeReader<FastaRecord> gnr({0,sg});
        FastaRecord node;
        while (gnr.next_record(node)) {
            kf.getMinSketch(node.seq, sketch);
            std::cout << "Node " << node.name << " sketch ";
            std::copy(sketch.cbegin(),sketch.cend(),std::ostream_iterator<MinPosIDX>(std::cout, "; "));
            std::cout << std::endl;
            for (const auto &sk : sketch) {
                kmer_to_graphposition[sk.hash].emplace_back(node.id, sk.pos);
            }
            sketch.clear();
        }
    }

    uint64_t map_reads(std::string filename, uint32_t error) {
        FastqReader<FastqRecord> fastqReader({0},filename);
        StrandedMinimiserSketchFactory kf(k, w);
        std::set<MinPosIDX> sketch;
        FastqRecord read;
        while(fastqReader.next_record(read)) {
            std::cout << "Mapping read " << read.name << std::endl;
            std::vector<Hit> hits;
            kf.getMinSketch(read.seq, sketch);
            for (const auto &sk:sketch){
                std::unordered_map<uint64_t, std::vector<graphStrandPos>>::const_iterator foundKey(kmer_to_graphposition.find(sk.hash));
                if (foundKey == kmer_to_graphposition.end()) continue;
                for (const auto &match : foundKey->second) {
                    if ((sk.pos ^ match.pos) >= 0) { // They are on the same strand (+ * -) or (- * +) if different strands
                        hits.emplace_back(match.node, true, sk.pos - match.pos, match.pos);
                    } else {
                        hits.emplace_back(match.node, false, sk.pos + match.pos, match.pos);
                    }
                }
            }

            std::sort(hits.begin(),hits.end());
            std::multiset<sgNodeID_t> nodesCt;
            for (const auto &h:hits) {
                nodesCt.insert(h.node);
            }
            for (const auto &nc: nodesCt){
                nodesCt.count(nc);
            }
            int b(0);
            for (unsigned int e = 0; e < hits.size(); e++) {
                if (e == hits.size() or
                        hits[e+1].node != hits[e].node or
                        hits[e+1].same_strand != hits[e].same_strand or
                        hits[e+1].intercept - hits[e].intercept >= error) {
                    // Insert maximal colineal subset hits[b,e] into C
                    std::copy(hits.begin()+b,hits.begin()+e, std::ostream_iterator<Hit>(std::cout,"; "));
                    std::cout << std::endl;
                    // Print left and right most query/target positions in C
                    b = e+1;
                }
            }
            sketch.clear();
        }
    }

    /*
    uint64_t map_reads(const uint16_t min_matches, std::string filename) {
        if (min_matches < 2 ) {
            throw std::invalid_argument("min_matches < 2, long reads can only be considered mapped if > 2 kmers are shared between read and sequence");
        }
        std::ofstream length_to_numMapped("ReadsNodesMapped.tsv");
        length_to_numMapped << "readLength\tnumMapped\n";
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
                mapping.read_id = read.id+1;
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

                    unsigned int mappedNodes(0);
                    for (auto &rm:read_mappings) {
                        if (rm.node != 0 and rm.unique_matches >= min_matches) {
#pragma omp critical(lr_reads_in_node)
                            {
                                if (nodes_in_read.size() < rm.read_id) { nodes_in_read.resize(rm.read_id + 100000); }
                                nodes_in_read[rm.read_id].push_back(rm.node);
                            }
                            mapped_count++;
                            mappedNodes++;
                        }

                    }
#pragma omp critical(log_mapped_count)
                    length_to_numMapped << read.seq.length() << "\t" << mappedNodes << "\n";
                }
                auto tc = ++total_count;
                if (tc % 100000 == 0) std::cout << mapped_count << " / " << tc << std::endl;
#pragma omp critical(lr_fastq_reader)
                {
                    c = fastqReader.next_record(read);
                }
            }

        }
        std::cout << int(k) << "-mers found " << kmers_found <<std::endl;
        std::cout<<"Reads mapped: "<<mapped_count<<" / "<<total_count<<std::endl;

        reads_in_node.resize(sg.nodes.size()+1);
        for (auto read = nodes_in_read.cbegin(); read != nodes_in_read.cend(); ++read){
            for (const auto &node:*read) {
                reads_in_node[std::abs(node)].push_back(read - nodes_in_read.cbegin());
            }
        }

        return mapped_count;
    }

*/
    uint64_t map_reads(std::unordered_set<uint64_t> const &  reads_to_remap={}) {

    }

    std::array<std::set<uint32_t>, 4> getReadSets(std::array<sgNodeID_t,4> &nodes) {
        std::array<std::set<uint32_t>, 4> readIDs;
        readIDs[0].insert(reads_in_node[std::abs(nodes[0])].cbegin(), reads_in_node[std::abs(nodes[0])].cend());
        readIDs[1].insert(reads_in_node[std::abs(nodes[1])].cbegin(), reads_in_node[std::abs(nodes[1])].cend());
        readIDs[2].insert(reads_in_node[std::abs(nodes[2])].cbegin(), reads_in_node[std::abs(nodes[2])].cend());
        readIDs[3].insert(reads_in_node[std::abs(nodes[3])].cbegin(), reads_in_node[std::abs(nodes[3])].cend());
        return readIDs;
    };
};


#endif //BSG_LONGREADMAPPER_H
