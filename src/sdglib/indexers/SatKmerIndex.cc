//
// Created by Luis Yanes (EI) on 2019-06-25.
//

#include "SatKmerIndex.hpp"
#include <sdglib/graph/SequenceDistanceGraph.hpp>

void SatKmerIndex::generate_index(const SequenceDistanceGraph &sg, uint8_t filter_limit, bool verbose) {
    // this can be parallelised by contig by aggregating different k_usage vectors on the first step, second and third steps are trickier
    // This two-pass approach could be easily adapted to a single memory block rather than vector of vectors
    uint64_t indexed_positions(0),filtered_kmers(0),total_kmers(0);
    //---- First Step, count k-mer occupancy ----//
    if (verbose) sdglib::OutputLog() << "First pass: Generating {kmer,contig,offset} " <<std::endl;
    kmerEnd.resize(std::pow(4,k));
    std::vector<kmerPos> all_kmers;
    all_kmers.reserve(100000000);
    std::vector<std::pair<bool,uint64_t > > contig_kmers;
    contig_kmers.reserve(1000000); //1Mbp per contig to start with?
    StringKMerFactory skf(k);
    for (sgNodeID_t n = 1; n < sg.nodes.size(); ++n) {
        if (sg.nodes[n].sequence.size() >= k) {
            contig_kmers.clear();
            skf.create_kmers(sg.nodes[n].sequence, contig_kmers);
            int k_i(0);
            for (const auto &kmer:contig_kmers) {
                all_kmers.emplace_back(kmer.second, n, kmer.first ? k_i + 1 : -(k_i + 1));
                k_i++;
            }
            total_kmers+=contig_kmers.size();
        }
    }
    //---- Second Step, reserve space for each vector in structure (avoiding reallocations and such)----//
    if (verbose) sdglib::OutputLog() << "Sorting, linearising structure and saving positions" << std::endl;
    sdglib::sort(all_kmers.begin(), all_kmers.end(), kmerPos::byKmerContigOffset());
    assembly_kmers.clear();

    if (verbose) sdglib::OutputLog() << "Filtering kmers appearing less than " << filter_limit << " from " << total_kmers << " initial kmers" << std::endl;
    auto witr = all_kmers.begin();
    auto ritr = witr;
    uint64_t ckmer(0);
    for (; ritr != all_kmers.end();) {

        while (ckmer != ritr->kmer) {
            ++ckmer;
            kmerEnd[ckmer]=kmerEnd[ckmer-1];
        }
        auto bitr = ritr;
        while (ritr != all_kmers.end() and bitr->kmer == ritr->kmer) {
            ++ritr;
        }
        if (ritr - bitr < filter_limit) {
            filtered_kmers+=1;
            kmerEnd[ckmer] += ritr-bitr;
            while (bitr != ritr) {
                contig_offsets.emplace_back(bitr->contigID, bitr->offset);
                ++witr;
                ++bitr;
            }
        }
    }
    while (ckmer < std::pow(4,k)) {
        ++ckmer;
        kmerEnd[ckmer] = kmerEnd[ckmer-1];
    }
    assert(kmerEnd[std::pow(4,k)-1] == contig_offsets.size());
    std::vector<kmerPos>().swap(all_kmers);
    if (verbose) sdglib::OutputLog() << "Kmers for mapping " << filtered_kmers << std::endl;
    if (verbose) sdglib::OutputLog() << "Elements in structure " << contig_offsets.size() << std::endl;
    if (verbose) sdglib::OutputLog() << "DONE" << std::endl;

}