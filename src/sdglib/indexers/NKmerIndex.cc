//
// Created by Luis Yanes (EI) on 2019-06-25.
//

#include "NKmerIndex.hpp"
#include <sdglib/graph/SequenceDistanceGraph.hpp>

uint64_t filter_kmers(std::vector<kmerPos> &kmers, int max_kmer_repeat) {
    uint64_t total_kmers(0);
    auto witr = kmers.begin();
    auto ritr = witr;
    for (; ritr != kmers.end();) {
        auto bitr = ritr;
        while (ritr != kmers.end() and bitr->kmer == ritr->kmer) {
            ++ritr;
        }
        if (ritr - bitr < max_kmer_repeat) {
            total_kmers+=1;
            while (bitr != ritr) {
                *witr = *bitr;
                ++witr;
                ++bitr;
            }
        }
    }
    kmers.resize(witr-kmers.begin());
    return total_kmers;
}

uint64_t filter_kmers(std::vector<kmerPos128> &kmers, int max_kmer_repeat) {
    uint64_t total_kmers(0);
    auto witr = kmers.begin();
    auto ritr = witr;
    for (; ritr != kmers.end();) {
        auto bitr = ritr;
        while (ritr != kmers.end() and bitr->kmer == ritr->kmer) {
            ++ritr;
        }
        if (ritr - bitr < max_kmer_repeat) {
            total_kmers+=1;
            while (bitr != ritr) {
                *witr = *bitr;
                ++witr;
                ++bitr;
            }
        }
    }
    kmers.resize(witr-kmers.begin());
    return total_kmers;
}

NKmerIndex::NKmerIndex(const SequenceDistanceGraph &_sg, uint8_t k, int filter_limit) : k(k), bfilter(70*1024*1024), sg(_sg){
    uint64_t total_length=0;
#pragma omp parallel for reduction(+:total_length)
    for (sgNodeID_t n = 1; n < sg.nodes.size(); ++n) {
        total_length+=sg.nodes[n].sequence.size();
    }
    assembly_kmers.reserve(total_length);
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
    auto total_kmers(filter_kmers(assembly_kmers, filter_limit));
#pragma omp parallel for
    for (uint64_t kidx = 0; kidx < assembly_kmers.size(); ++kidx) {
        bfilter.add(assembly_kmers[kidx].kmer);
    }
}

void NKmerIndex::get_all_kmer_matches(std::vector<std::vector<std::pair<int32_t, int32_t>>> &matches,
                                      std::vector<std::pair<bool, uint64_t>> &seq_kmers) {
    if (matches.size() < seq_kmers.size()) matches.resize(seq_kmers.size());
    for (auto i = 0; i < seq_kmers.size(); ++i) {
        matches[i].clear();
        auto first = find(seq_kmers[i].second);
        for (auto it = first; it != assembly_kmers.end() && it->kmer == seq_kmers[i].second; ++it) {
            int32_t offset = it->offset; //so far, this is +1 and the sign indicates direction of kmer in contig
            sgNodeID_t node = it->contigID; //so far, this is always positive
            if (seq_kmers[i].first != (offset > 0)) { //match is on reverse
                node = -node;
                offset = ((int) sg.nodes[std::llabs(it->contigID)].sequence.size()) - std::abs(offset);
            } else offset = std::abs(offset) - 1;
            matches[i].emplace_back(node, offset);
        }
    }
}

NKmerIndex128::NKmerIndex128(const SequenceDistanceGraph &_sg, uint8_t k, int filter_limit) : k(k), bfilter(70*1024*1024), sg(_sg){
    uint64_t total_length=0;
#pragma omp parallel for reduction(+:total_length)
    for (sgNodeID_t n = 1; n < sg.nodes.size(); ++n) {
        total_length+=sg.nodes[n].sequence.size();
    }
    assembly_kmers.reserve(total_length);
#pragma omp parallel
    {
        StringKMerFactory128 skf(k);
        std::vector<kmerPos128> local_kmers;
        local_kmers.reserve(total_length/(omp_get_num_threads()*4));
        std::vector<std::pair<bool,__uint128_t > > contig_kmers;
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


    sdglib::sort(assembly_kmers.begin(),assembly_kmers.end(), kmerPos128::byKmerContigOffset());
    auto total_kmers(filter_kmers(assembly_kmers, filter_limit));
#pragma omp parallel for
    for (uint64_t kidx = 0; kidx < assembly_kmers.size(); ++kidx) {
        bfilter.add(assembly_kmers[kidx].kmer);
    }
}

void NKmerIndex128::get_all_kmer_matches(std::vector<std::vector<std::pair<int32_t, int32_t>>> &matches,
                                      std::vector<std::pair<bool, __uint128_t>> &seq_kmers) {
    if (matches.size() < seq_kmers.size()) matches.resize(seq_kmers.size());
    for (auto i = 0; i < seq_kmers.size(); ++i) {
        matches[i].clear();
        auto first = find(seq_kmers[i].second);
        for (auto it = first; it != assembly_kmers.end() && it->kmer == seq_kmers[i].second; ++it) {
            int32_t offset = it->offset; //so far, this is +1 and the sign indicates direction of kmer in contig
            sgNodeID_t node = it->contigID; //so far, this is always positive
            if (seq_kmers[i].first != (offset > 0)) { //match is on reverse
                node = -node;
                offset = ((int) sg.nodes[std::llabs(it->contigID)].sequence.size()) - std::abs(offset);
            } else offset = std::abs(offset) - 1;
            matches[i].emplace_back(node, offset);
        }
    }
}