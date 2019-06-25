//
// Created by Luis Yanes (EI) on 2019-06-25.
//

#include "SatKmerIndex.hpp"
#include <sdglib/graph/SequenceDistanceGraph.hpp>

void SatKmerIndex::generate_index_vec_allocs(const SequenceDistanceGraph &sg, int filter_limit, bool verbose) {
    uint64_t total_kmers(0);
    assembly_kmers.resize(std::pow(4,k));
    if (verbose) {
        sdglib::OutputLog() << "Updating mapping index for k=" << std::to_string(k) << std::endl;
        sdglib::OutputLog() << "Number of kmers to store " << std::to_string(std::pow(4,k)) << std::endl;
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

    if (verbose) sdglib::OutputLog() << "Filtering kmers appearing less than " << filter_limit << " from " << total_kmers << " initial kmers" << std::endl;
    auto pair = filter_kmers(filter_limit);
    if (verbose) sdglib::OutputLog() << "Kmers for mapping " << pair.first << std::endl;
    if (verbose) sdglib::OutputLog() << "Elements in structure " << pair.second << std::endl;
    if (verbose) sdglib::OutputLog() << "DONE" << std::endl;
}

void SatKmerIndex::generate_index_prealloc(const SequenceDistanceGraph &sg, uint8_t filter_limit, bool verbose) {
    // this can be parallelised by contig by aggregating different k_usage vectors on the first step, second and third steps are trickier
    // This two-pass approach could be easily adapted to a single memory block rather than vector of vectors
    uint64_t indexed_positions(0),indexed_kmers(0),total_kmers(0);
    //---- First Step, count k-mer occupancy ----//
    if (verbose) sdglib::OutputLog() << "First pass: analysing usage over " << std::to_string(uint64_t (std::pow(4,k))) << " k-mers..."<<std::endl;
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
    if (verbose) sdglib::OutputLog() << "Reserving space on internal vectors..." << std::endl;
    assembly_kmers.clear();
    assembly_kmers.reserve(k_usage.size());
    for (auto &count:k_usage){
        if (count<filter_limit) { //Already accounts for the filter
            assembly_kmers.emplace_back(count);
            ++indexed_kmers;
        }
        else{
            assembly_kmers.emplace_back();
        }
    }
    if (verbose) sdglib::OutputLog() << "Space reserved to index "<< indexed_kmers <<" k-mers with less than "<<filter_limit+1<<" copies..." << std::endl;

    //---- Third step, actually save the positions ---//
    if (verbose) sdglib::OutputLog() << "Second pass: filling index..." << std::endl;
    contig_kmers.reserve(10000);
    for (sgNodeID_t n = 1; n < sg.nodes.size(); ++n) {
        if (sg.nodes[n].sequence.size() >= k) {
            contig_kmers.clear();
            skf.create_kmers(sg.nodes[n].sequence, contig_kmers);
            int k_i(0);
            for (const auto &kmer:contig_kmers) {
                if (k_usage[kmer.second]<filter_limit){ //filter, only inserts on relevant k-mers, which are already reserved too
                    assembly_kmers[kmer.second].emplace_back(n, kmer.first ? k_i + 1 : -(k_i + 1)); //is this really how we save the offsets on the indexes?
                    ++indexed_positions;
                }
                ++k_i;
            }
            total_kmers+=contig_kmers.size();
        }
    }

    if (verbose) sdglib::OutputLog() << indexed_positions << "/" << total_kmers << " postions indexed by " << indexed_kmers << " kmers" << std::endl;

}

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

void SatKmerIndex::generate_index_parallel(const SequenceDistanceGraph &sg, int filter_limit, bool verbose) {
    uint64_t total_kmers(0);
    assembly_kmers.resize(std::pow(4,k));
    if (verbose) {
        sdglib::OutputLog() << "Updating mapping index for k=" << std::to_string(k) << std::endl;
        sdglib::OutputLog() << "Number of kmers to store " << std::to_string(uint64_t (std::pow(4,k))) << std::endl;
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

    if (verbose) sdglib::OutputLog() << "Filtering kmers appearing less than " << filter_limit << " from " << total_kmers << " initial kmers" << std::endl;
    auto pair = filter_kmers(filter_limit);
    if (verbose) sdglib::OutputLog() << "Kmers for mapping " << pair.first << std::endl;
    if (verbose) sdglib::OutputLog() << "Elements in structure " << pair.second << std::endl;
    if (verbose) sdglib::OutputLog() << "DONE" << std::endl;
}
