//
// Created by Luis Yanes (EI) on 17/08/2018.
//

#ifndef BSG_SATKMERINDEX_HPP
#define BSG_SATKMERINDEX_HPP

#include <vector>
#include <sglib/utilities/omp_safe.hpp>
#include <sglib/types/KmerTypes.hpp>
#include <sglib/factories/KMerFactory.hpp>
#include <sglib/logger/OutputLog.hpp>
#include <sglib/graph/SequenceGraph.hpp>

struct ContigOffset {
    ContigOffset() = default;
    ContigOffset(uint32_t contigID, int32_t offset) : contigID(contigID),offset(offset) {}

    int32_t contigID = 0;
    int32_t offset = 0;

    const bool operator==(const kmerPos &a) const { return std::tie(contigID, offset) == std::tie(a.contigID, a.offset);}
};

class SatKmerIndex {
    std::vector<std::vector<ContigOffset>> assembly_kmers;
    uint8_t k=15;
public:
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

    void generate_index_vec_allocs(const SequenceGraph &sg, int filter_limit = 200, bool verbose=true) {
        uint64_t total_kmers(0);
        assembly_kmers.resize(std::pow(4,k));
        if (verbose) {
            sglib::OutputLog() << "Updating mapping index for k=" << std::to_string(k) << std::endl;
            sglib::OutputLog() << "Number of kmers to store " << std::to_string(std::pow(4,k)) << std::endl;
        }

        StringKMerFactory skf(k);
        std::vector<std::pair<bool,uint64_t > > contig_kmers;
        contig_kmers.reserve(10000);
        for (sgNodeID_t n = 1; n < sg.nodes.size(); ++n) {
            if (sg.nodes[n].sequence.size() >= k) {
                contig_kmers.clear();
                skf.create_kmers(sg.nodes[n].sequence, contig_kmers);
                int k_i(0);
                for (const auto &kmer:contig_kmers) {
                    if (assembly_kmers[kmer.second].size()<filter_limit+1) assembly_kmers[kmer.second].emplace_back(n, kmer.first ? k_i + 1 : -(k_i + 1));
                    k_i++;
                }
                total_kmers+=contig_kmers.size();
            }
        }

        if (verbose) sglib::OutputLog() << "Filtering kmers appearing less than " << filter_limit << " from " << total_kmers << " initial kmers" << std::endl;
        auto pair = filter_kmers(filter_limit);
        if (verbose) sglib::OutputLog() << "Kmers for mapping " << pair.first << std::endl;
        if (verbose) sglib::OutputLog() << "Elements in structure " << pair.second << std::endl;
        if (verbose) sglib::OutputLog() << "DONE" << std::endl;
    }

    /**
     * This function is a alternative generator of the kmer index, reserving space first in the internal vectors.
     * @param sg
     * @param filter_limit
     * @param verbose
     */
    void generate_index(const SequenceGraph &sg, uint8_t filter_limit = 200, bool verbose=true) {
        // this can be parallelised by contig by aggregating different k_usage vectors on the first step, second and third steps are trickier
        // This two-pass approach could be easily adapted to a single memory block rather than vector of vectors
        uint64_t indexed_positions(0),indexed_kmers(0),total_kmers(0);
        //---- First Step, count k-mer occupancy ----//
        if (verbose) sglib::OutputLog() << "First pass: analysing usage over " << std::to_string(uint64_t (std::pow(4,k))) << " k-mers..."<<std::endl;
        std::vector<uint8_t> k_usage(std::pow(4,k));
        std::vector<std::pair<bool,uint64_t > > contig_kmers;
        contig_kmers.reserve(1000000); //1Mbp per contig to start with?
        StringKMerFactory skf(k);
        for (sgNodeID_t n = 1; n < sg.nodes.size(); ++n) {
            if (sg.nodes[n].sequence.size() >= k) {
                contig_kmers.clear();
                skf.create_kmers(sg.nodes[n].sequence, contig_kmers);
                for (const auto &kmer:contig_kmers) if (k_usage[kmer.second]<=filter_limit) ++k_usage[kmer.second]; //beware of overflow, hence stop counting after reaching limit
            }
        }
        //---- Second Step, reserve space for each vector in structure (avoiding reallocations and such)----//
        if (verbose) sglib::OutputLog() << "Reserving space on internal vectors..." << std::endl;
        assembly_kmers.clear();
        assembly_kmers.reserve(k_usage.size());
        for (auto &count:k_usage){
            if (count<=filter_limit) { //Already accounts for the filter
                assembly_kmers.emplace_back(count);
                ++indexed_kmers;
            }
            else{
                assembly_kmers.emplace_back();
            }
        }
        if (verbose) sglib::OutputLog() << "Space reserved to index "<< indexed_kmers <<" k-mers with less than "<<filter_limit+1<<" copies..." << std::endl;

        //---- Third step, actually save the positions ---//
        if (verbose) sglib::OutputLog() << "Second pass: filling index..." << std::endl;
        contig_kmers.reserve(10000);
        for (sgNodeID_t n = 1; n < sg.nodes.size(); ++n) {
            if (sg.nodes[n].sequence.size() >= k) {
                contig_kmers.clear();
                skf.create_kmers(sg.nodes[n].sequence, contig_kmers);
                int k_i(0);
                for (const auto &kmer:contig_kmers) {
                    if (k_usage[kmer.second]<=filter_limit){ //filter, only inserts on relevant k-mers, which are already reserved too
                        assembly_kmers[kmer.second].emplace_back(n, kmer.first ? k_i + 1 : -(k_i + 1)); //is this really how we save the offsets on the indexes?
                        ++indexed_positions;
                    }
                    ++k_i;
                }
                total_kmers+=contig_kmers.size();
            }
        }

        if (verbose) sglib::OutputLog() << indexed_positions << "/" << total_kmers << " postions indexed by " << indexed_kmers << " kmers" << std::endl;

    }

    void generate_index_parallel(const SequenceGraph &sg, int filter_limit = 200, bool verbose=true) {
        uint64_t total_kmers(0);
        assembly_kmers.resize(std::pow(4,k));
        if (verbose) {
            sglib::OutputLog() << "Updating mapping index for k=" << std::to_string(k) << std::endl;
            sglib::OutputLog() << "Number of kmers to store " << std::to_string(uint64_t (std::pow(4,k))) << std::endl;
        }
#pragma omp parallel reduction(+:total_kmers)
        {
            StringKMerFactory skf(k);
            std::vector<kmerPos> local_kmers;
            local_kmers.reserve(80000);
            std::vector<std::pair<bool,uint64_t > > contig_kmers;
            contig_kmers.reserve(10000);
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
                    total_kmers+=contig_kmers.size();
                }

                if (local_kmers.size() > 8000000) {
#pragma omp critical(push_kmers)
                    {
                        for (const auto &lk : local_kmers) {
                            assembly_kmers[lk.kmer].emplace_back(lk.contigID,lk.offset);
                        }
                    }
                    local_kmers.clear();
                }
            }
            if (!local_kmers.empty()) {
#pragma omp critical(push_kmers)
                {
                    for (const auto &lk : local_kmers) {
                        assembly_kmers[lk.kmer].emplace_back(lk.contigID,lk.offset);
                    }
                }
                local_kmers.clear();
            }
        }

        if (verbose) sglib::OutputLog() << "Filtering kmers appearing less than " << filter_limit << " from " << total_kmers << " initial kmers" << std::endl;
        auto pair = filter_kmers(filter_limit);
        if (verbose) sglib::OutputLog() << "Kmers for mapping " << pair.first << std::endl;
        if (verbose) sglib::OutputLog() << "Elements in structure " << pair.second << std::endl;
        if (verbose) sglib::OutputLog() << "DONE" << std::endl;
    }

    bool empty(uint64_t kmer) const { return assembly_kmers[kmer].empty(); }
    const_iterator begin(uint64_t kmer) const {return assembly_kmers[kmer].cbegin();}
    const_iterator end(uint64_t kmer) const {return assembly_kmers[kmer].cend();}

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
