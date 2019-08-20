//
// Created by Luis Yanes (EI) on 2019-07-08.
//

#include <catch.hpp>
#include <unordered_map>
#include <sdglib/batch_counter/BatchKmersCounter.hpp>
#include <sdglib/workspace/WorkSpace.hpp>
#include <cstring>

TEST_CASE("Counts 60-mers correctly from PairedReadDatastore") {
    {
        std::string r1_filepath("../tests/datasets/workspace/pe/pe_R1.fastq");
        std::string r2_filepath("../tests/datasets/workspace/pe/pe_R2.fastq");
        std::string prds_output_path("pe.prseq");
        PairedReadsDatastore::build_from_fastq(prds_output_path, r1_filepath, r2_filepath, prds_output_path, 0, 350);
    }

    WorkSpace ws;
    const PairedReadsDatastore ds(ws,"pe.prseq");

    auto kmer_list = BatchKmersCounter::buildKMerCount(60, ds, 0, "./", "./", 0);

    __uint128_t kmer;
    std::vector<std::pair<__uint128_t, uint64_t>> batch_count;
    batch_count.resize(kmer_list->size);
    for (uint64_t i =0; i < kmer_list->size; i++) {
        memcpy(&kmer, &kmer_list->kmers[i], 16);
        batch_count[i] = std::make_pair(kmer, kmer_list->kmers[i].count);
    }

    sdglib::sort(batch_count.begin(), batch_count.end(), [](const std::pair<__uint128_t, uint64_t> &a, const std::pair<__uint128_t, uint64_t> &b){return a.first < b.first;});

    std::string filepath("../tests/datasets/workspace/pe/60mer_counts");
    std::ifstream counts_file(filepath);
    if (!counts_file) {
        std::cerr << "Error opening " << filepath << ": " << std::strerror(errno) << std::endl;
        throw std::runtime_error("Error opening " + filepath + ": " + std::strerror(errno));
    }

    StreamKmerFactory128 skf(60);
    std::string skmer;
    std::string count;
    std::vector<std::pair<__uint128_t, uint64_t>> file_count;
    file_count.reserve(kmer_list->size);
    std::vector<__uint128_t> kmers;
    while (counts_file >> skmer >> count) {
        skf.produce_all_kmers(skmer.c_str(), kmers);
        file_count.emplace_back(kmers[0], std::min(uint64_t(255), std::stoull(count.data())));
        kmers.clear();
    }

    REQUIRE(batch_count == file_count);

    ::unlink("small_K.freqs");
    ::unlink("pe.prseq");

}

TEST_CASE("Counts 60-mers correctly from PairedReadDatastore, using batches") {
    {
        std::string r1_filepath("../tests/datasets/workspace/pe/pe_R1.fastq");
        std::string r2_filepath("../tests/datasets/workspace/pe/pe_R2.fastq");
        std::string prds_output_path("pe.prseq");
        PairedReadsDatastore::build_from_fastq(prds_output_path, r1_filepath, r2_filepath, prds_output_path, 0, 350);
    }

    WorkSpace ws;
    const PairedReadsDatastore ds(ws,"pe.prseq");

    auto kmer_list = BatchKmersCounter::buildKMerCount(60, ds, 0, "./", "./", 4);

    std::vector<std::pair<__uint128_t, uint64_t>> batch_count;
    __uint128_t kmer;
    batch_count.resize(kmer_list->size);
    for (uint64_t i =0; i < kmer_list->size; i++) {
        memcpy(&kmer, &kmer_list->kmers[i], 16);
        batch_count[i] = std::make_pair(kmer, kmer_list->kmers[i].count);
    }

    sdglib::sort(batch_count.begin(), batch_count.end(), [](const std::pair<__uint128_t, uint64_t> &a, const std::pair<__uint128_t, uint64_t> &b){return a.first < b.first;});


    std::string filepath("../tests/datasets/workspace/pe/60mer_counts");
    std::ifstream counts_file(filepath);
    if (!counts_file) {
        std::cout << "Error opening " << filepath << ": " << std::strerror(errno) << std::endl;
        throw std::runtime_error("Error opening " + filepath + ": " + std::strerror(errno));
    }

    StreamKmerFactory128 skf(60);
    std::string skmer;
    std::string count;
    std::vector<std::pair<__uint128_t, uint64_t>> file_count;
    file_count.reserve(kmer_list->size);
    std::vector<__uint128_t> kmers;
    while (counts_file >> skmer >> count) {
        skf.produce_all_kmers(skmer.c_str(), kmers);
        file_count.emplace_back(kmers[0], std::min(255ul, std::stoul(count.data())));
        kmers.clear();
    }

    REQUIRE( batch_count.size() == file_count.size() );
    REQUIRE(batch_count == file_count);
    ::unlink("small_K.freqs");
    ::unlink("pe.prseq");

}