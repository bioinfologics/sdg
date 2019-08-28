//
// Created by Luis Yanes (EI) on 17/08/2018.
//

#pragma once

#include <cassert>
#include <vector>
#include <sdglib/utilities/omp_safe.hpp>
#include <sdglib/types/KmerTypes.hpp>
#include <sdglib/factories/KMerFactory.hpp>

class SequenceDistanceGraph;
struct ContigOffset {
    ContigOffset() = default;
    ContigOffset(uint32_t contigID, int32_t offset, int32_t rcOffset) : contigID(contigID),offset(offset),rcOffset(rcOffset) {}

    int32_t contigID = 0;
    int32_t offset = 0;
    int32_t rcOffset = 0;

    const bool operator==(const kmerPos &a) const { return std::tie(contigID, offset) == std::tie(a.contigID, a.offset);}
};

class SatKmerIndex {
    std::vector<uint64_t> kmerEnd;
    uint8_t k=0;
public:
    std::vector<ContigOffset> contig_offsets;
    using const_iterator = std::vector<ContigOffset>::const_iterator;

    SatKmerIndex(){};
    SatKmerIndex(const SequenceDistanceGraph &sg, uint8_t k=15, uint8_t filter_limit = 200);

    uint64_t beginCO(uint64_t kmer) const { return (0ull==kmer) ? 0ull : kmerEnd[kmer-1];}
    uint64_t endCO(uint64_t kmer) const { return kmerEnd[kmer]; }
};
