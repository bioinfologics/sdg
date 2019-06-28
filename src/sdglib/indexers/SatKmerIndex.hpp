//
// Created by Luis Yanes (EI) on 17/08/2018.
//

#ifndef BSG_SATKMERINDEX_HPP
#define BSG_SATKMERINDEX_HPP

#include <cassert>
#include <vector>
#include <sdglib/utilities/omp_safe.hpp>
#include <sdglib/types/KmerTypes.hpp>
#include <sdglib/factories/KMerFactory.hpp>
#include <sdglib/logger/OutputLog.hpp>

class SequenceDistanceGraph;
struct ContigOffset {
    ContigOffset() = default;
    ContigOffset(uint32_t contigID, int32_t offset) : contigID(contigID),offset(offset) {}

    int32_t contigID = 0;
    int32_t offset = 0;

    const bool operator==(const kmerPos &a) const { return std::tie(contigID, offset) == std::tie(a.contigID, a.offset);}
};

class SatKmerIndex {
    std::vector<std::vector<ContigOffset>> assembly_kmers;
    std::vector<uint64_t> kmerEnd;
    uint8_t k=15;
public:
    std::vector<ContigOffset> contig_offsets;
    using const_iterator = std::vector<ContigOffset>::const_iterator;

    SatKmerIndex(){}
    explicit SatKmerIndex(uint8_t k) : k(k) {
        if (k > 15) {
            throw std::runtime_error(
                    "You are trying to use K>15, which is not supported by this structure. "
                    "Please consider NKmerIndex or UniqueKmerIndex as alternatives");
        }
    }

    std::pair<uint64_t, uint64_t> filter_kmers(int max_kmer_repeat) {
        uint64_t num_kmers=0, num_elements=0;
#pragma omp parallel for reduction(+:num_kmers, num_elements)
        for (uint64_t kidx=0; kidx < assembly_kmers.size(); ++kidx) {
            if (assembly_kmers[kidx].size() >= max_kmer_repeat) {
                std::vector<ContigOffset>().swap(assembly_kmers[kidx]);
            }
            if (!assembly_kmers[kidx].empty()){num_kmers++; num_elements+=assembly_kmers[kidx].size();}
        }

        return {num_kmers, num_elements};
    }

    /**
     * @brief
     * Generate an index for the start location of each kmer and a list of {+-node,pos} that can be queried by kmer.
     * @param sg
     * @param filter_limit
     * @param verbose
     */
    void generate_index(const SequenceDistanceGraph &sg, uint8_t filter_limit = 200, bool verbose=true);

    bool empty(uint64_t kmer) const { return assembly_kmers[kmer].empty(); }
    const_iterator begin(uint64_t kmer) const {return assembly_kmers[kmer].cbegin();}
    const_iterator end(uint64_t kmer) const {return assembly_kmers[kmer].cend();}

    auto beginCO(uint64_t kmer) const { return (0ull==kmer) ? 0ull : kmerEnd[kmer-1];}
    auto endCO(uint64_t kmer) const { return kmerEnd[kmer]; }

    /**
     * @brief
     * Returns the index to the first and last [+-node,pos] in this kmer
     * @param kmer Query kmer
     * @return
     */
    const std::vector<ContigOffset> & find(const uint64_t kmer) const {
        return assembly_kmers[kmer];
    }
};


#endif //BSG_SATKMERINDEX_HPP
