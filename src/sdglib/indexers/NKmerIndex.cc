//
// Created by Luis Yanes (EI) on 2019-06-25.
//

#include "NKmerIndex.hpp"
#include <sdglib/graph/SequenceDistanceGraph.hpp>

void NKmerIndex::generate_index(const SequenceDistanceGraph &sg, int filter_limit, bool verbose) {
    uint64_t total_length=0;
#pragma omp parallel for reduction(+:total_length)
    for (sgNodeID_t n = 1; n < sg.nodes.size(); ++n) {
        total_length+=sg.nodes[n].sequence.size();
    }
    assembly_kmers.reserve(total_length);
    if (verbose) sdglib::OutputLog() << "Updating mapping index for k=" << std::to_string(k) << std::endl;
#pragma omp parallel
    {
        StringKMerFactory skf(k);
        std::vector<kmerPos> local_kmers;
        local_kmers.reserve(total_length/(omp_get_num_threads()*4));
        std::vector<std::pair<bool,uint64_t > > contig_kmers;
        contig_kmers.reserve(10000000);
#pragma omp for
        for (sgNodeID_t n = 1; n < sg.nodes.size(); ++n) {
            if (sg.nodes[n].sequence.size() >= k) {
                contig_kmers.clear();
                skf.create_kmers(sg.nodes[n].sequence, contig_kmers);
                int k_i(0);
                for (const auto &kmer:contig_kmers) {
                    local_kmers.emplace_back(kmer.second, n, kmer.first ? k_i + 1 : -(k_i + 1));
                    k_i++;
                }
            }

            if (local_kmers.size() > 800000000) {
#pragma omp critical(push_kmers)
                {
                    assembly_kmers.insert(assembly_kmers.end(), local_kmers.begin(), local_kmers.end());
                }
                local_kmers.clear();
            }
        }
        if (!local_kmers.empty()) {
#pragma omp critical(push_kmers)
            {
                assembly_kmers.insert(assembly_kmers.end(), local_kmers.begin(), local_kmers.end());
            }
            local_kmers.clear();
        }
    }


    sdglib::sort(assembly_kmers.begin(),assembly_kmers.end(), kmerPos::byKmerContigOffset());

    if (verbose) sdglib::OutputLog() << "Filtering kmers appearing less than " << filter_limit << " from " << assembly_kmers.size() << " initial kmers" << std::endl;
    auto total_kmers(filter_kmers(assembly_kmers, filter_limit));
#pragma omp parallel for
    for (uint64_t kidx = 0; kidx < assembly_kmers.size(); ++kidx) {
        filter.add(assembly_kmers[kidx].kmer);
    }

    if (verbose) sdglib::OutputLog() << "Elements for mapping " << assembly_kmers.size() << std::endl;
    if (verbose) sdglib::OutputLog() << "Total distinct kmers " << total_kmers << std::endl;
    if (verbose) sdglib::OutputLog() << "Number of elements in bloom " << filter.number_bits_set() << std::endl;
    if (verbose) sdglib::OutputLog() << "Filter FPR " << filter.false_positive_rate() << std::endl;
    if (verbose) sdglib::OutputLog() << "DONE" << std::endl;
}
