//
// Created by Luis Yanes (EI) on 17/08/2018.
//

#ifndef BSG_NKMERINDEX_HPP
#define BSG_NKMERINDEX_HPP

#include <vector>
#include <sglib/utilities/omp_safe.hpp>
#include <sglib/types/KmerTypes.hpp>
#include <sglib/factories/KMerFactory.hpp>
#include <sglib/logger/OutputLog.hpp>
#include <sglib/graph/SequenceGraph.hpp>
#include <sglib/bloom/BloomFilter.hpp>

class NKmerIndex {
    BloomFilter filter;
    std::vector<kmerPos> assembly_kmers;
    uint8_t k;
public:
    using const_iterator = std::vector<kmerPos>::const_iterator;

    explicit NKmerIndex(uint8_t k) : k(k), filter(70*1024*1024) {}

    void generate_index(const SequenceGraph &sg, int filter_limit = 200, bool verbose=true) {
        assembly_kmers.reserve(100000000);

        if (verbose) sglib::OutputLog() << "Updating mapping index for k=" << std::to_string(k) << std::endl;
#pragma omp parallel
        {
            StringKMerFactory skf(k);
            std::vector<kmerPos> local_kmers;
            local_kmers.reserve(10000000);
            std::vector<std::pair<bool,uint64_t > > contig_kmers;
            contig_kmers.reserve(1000000);
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

                if (local_kmers.size() > 8000000) {
#pragma omp critical(push_kmers)
                    {
                        assembly_kmers.insert(assembly_kmers.end(), local_kmers.begin(), local_kmers.end());
                        local_kmers.clear();
                    }
                }
            }
            if (!local_kmers.empty()) {
#pragma omp critical(push_kmers)
                {
                    assembly_kmers.insert(assembly_kmers.end(), local_kmers.begin(), local_kmers.end());
                    local_kmers.clear();
                }
            }
        }


#ifdef _OPENMP
        __gnu_parallel::sort(assembly_kmers.begin(),assembly_kmers.end(), kmerPos::byKmerContigOffset());
#else
        std::sort(assembly_kmers.begin(),assembly_kmers.end(), kmerPos::byKmerContigOffset());
#endif

        int max_kmer_repeat(filter_limit);
        if (verbose) sglib::OutputLog() << "Filtering kmers appearing less than " << max_kmer_repeat << " from " << assembly_kmers.size() << " initial kmers" << std::endl;
        auto witr = assembly_kmers.begin();
        auto ritr = witr;
        for (; ritr != assembly_kmers.end();) {
            auto bitr = ritr;
            while (ritr != assembly_kmers.end() and bitr->kmer == ritr->kmer) {
                ++ritr;
            }
            if (ritr-bitr < max_kmer_repeat) {
                while (bitr != ritr) {
                    filter.add(bitr->kmer);
                    *witr = *bitr;
                    ++witr;++bitr;
                }
            }
        }
        assembly_kmers.resize(witr-assembly_kmers.begin());

        if (verbose) sglib::OutputLog() << "Kmers for mapping " << assembly_kmers.size() << std::endl;
        if (verbose) sglib::OutputLog() << "Number of elements in bloom " << filter.number_bits_set() << std::endl;
        if (verbose) sglib::OutputLog() << "Filter FPR " << filter.false_positive_rate() << std::endl;
        if (verbose) sglib::OutputLog() << "DONE" << std::endl;
    }

    bool empty() const { return assembly_kmers.empty(); }
    const_iterator begin() const {return assembly_kmers.cbegin();}
    const_iterator end() const {return assembly_kmers.cend();}

    const_iterator find(const uint64_t kmer) const {
        if (filter.contains(kmer)) { return std::lower_bound(assembly_kmers.cbegin(), assembly_kmers.cend(), kmer); }
        return assembly_kmers.cend();
//        return std::lower_bound(assembly_kmers.cbegin(), assembly_kmers.cend(), kmer);
    }
};


#endif //BSG_NKMERINDEX_HPP
