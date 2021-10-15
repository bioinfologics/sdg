//
// Created by Luis Yanes (EI) on 17/08/2018.
//

#pragma once

#include <vector>
#include <sdglib/utilities/omp_safe.hpp>
#include <sdglib/types/KmerTypes.hpp>
#include <sdglib/factories/KMerFactory.hpp>
#include <sdglib/utilities/OutputLog.hpp>
#include <sdglib/graph/SequenceDistanceGraph.hpp>
#include <sdglib/bloom/BloomFilter.hpp>

class SequenceDistanceGraph;

/** @brief Indexes the graph by kmers allowing repeated kmers and provides the functions to operate with the index.
 *
 * Usually employed at the start of a mapping process.
 *
 * The class keeps 2 collections the assembly_kmers (vector<kmerPos>) storing kmer information and a
 * bfilter (BloomFilter) for quick kmer presence/absence checking and high frequency filtering.
 *
 */
class NKmerIndex {
    BloomFilter bfilter;
    std::vector<kmerPos> assembly_kmers;
    uint8_t k=0;
    const SequenceDistanceGraph &sg;
public:
    using const_iterator = std::vector<kmerPos>::const_iterator;

    /** @brief Initializes the index
     * @param _sg graph to index
     * @param k kmer size to use in the index
     * @param filter_limit maximum frequency of a kmer to be included in the graph ("bloom filter" threshold)
     */
    explicit NKmerIndex(const SequenceDistanceGraph &_sg,uint8_t k=15, int filter_limit = 200);

    /** @brief Initializes the index
     * @param _sg graph to index
     * @param k kmer size to use in the index
     * @param whitelist a list of kmers to restrict the index
     */
    explicit NKmerIndex(const SequenceDistanceGraph &_sg,uint8_t k, std::unordered_set<uint64_t> whitelist);

    /** @Empties the filter (only the assembly_kmers vector)
     *  TODO: shouldn't this empty the bfilter as well!?
     *
     * @return
     */
    bool empty() const { return assembly_kmers.empty(); }

    const_iterator begin() const {return assembly_kmers.cbegin();}
    const_iterator end() const {return assembly_kmers.cend();}

    /** @brief check for a kmer in the filter
     *
     * @param kmer int representation of a kmer
     * @return true if the kmer is in the filter false otherwise
     */
    bool filter(const uint64_t kmer) const {
        return bfilter.contains(kmer);
    }

    /** @brief find a kmer in the filter
     *
     * @param kmer int representation of a kmer
     * @return if the kmer was found an iterator to the first ocurrence of the kmer otherwise .end() (as in .find())
     */
    const_iterator find(const uint64_t kmer) const {
        if (bfilter.contains(kmer)) { return std::lower_bound(assembly_kmers.cbegin(), assembly_kmers.cend(), kmer); }
        return assembly_kmers.cend();
    }

    /** @brief Get all matches to the sequence kmers queried
     *
     * Given a vector of kmers coming from the sequence returns all matches to those kmer in the index.
     * TODO: shouldn't this return a kmerpos like struct like the other functions that deal with kmer positions instead of a node id and a node offset
     * TODO: sus int32_t tratment of node ids :/
     *
     * @param matches empty vector to put the result in
     * @param seq_kmers sequence kmers to be queried in the index
     */
    void get_all_kmer_matches(std::vector<std::vector<std::pair<int32_t, int32_t>>> & matches, std::vector<std::pair<bool, uint64_t>> & seq_kmers);

};


class NKmerIndex128 {
    BloomFilter bfilter;
    std::vector<kmerPos128> assembly_kmers;
    uint8_t k=0;
    const SequenceDistanceGraph &sg;
public:
    using const_iterator = std::vector<kmerPos128>::const_iterator;

    explicit NKmerIndex128(const SequenceDistanceGraph &_sg,uint8_t k=63, int filter_limit = 200);

    bool empty() const { return assembly_kmers.empty(); }
    const_iterator begin() const {return assembly_kmers.cbegin();}
    const_iterator end() const {return assembly_kmers.cend();}

    bool filter(const __uint128_t kmer) const {
        return bfilter.contains(kmer);
    }

    const_iterator find(const __uint128_t kmer) const {
        if (bfilter.contains(kmer)) { return std::lower_bound(assembly_kmers.cbegin(), assembly_kmers.cend(), kmer); }
        return assembly_kmers.cend();
    }

    void get_all_kmer_matches(std::vector<std::vector<std::pair<int32_t, int32_t>>> & matches, std::vector<std::pair<bool, __uint128_t>> & seq_kmers);

};
