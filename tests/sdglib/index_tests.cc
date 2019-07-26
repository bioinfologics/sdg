//
// Created by Luis Yanes (EI) on 09/10/2018.
//

#include <catch.hpp>
#include <sdglib/indexers/UniqueKmerIndex.hpp>
#include <sdglib/indexers/NKmerIndex.hpp>
#include <sdglib/workspace/WorkSpace.hpp>

TEST_CASE("UniqueKmerIndex create and lookup") {
    unsigned int K(15);
    std::vector<KmerIDX> readkmers;
    StreamKmerIDXFactory skf(K);
    std::string seqMissing = "AAAAAAAAAAAAAAA";
    std::string seqPresent = "CTTGCGGGTTTCCAG";
    WorkSpace ws;
    SequenceDistanceGraph sg(ws);
    sg.load_from_gfa("../tests/datasets/graph/tgraph.gfa");
    UniqueKmerIndex ukm(sg, K);

    REQUIRE(ukm.getMap().size() != 0);

    skf.produce_all_kmers(seqMissing.data(),readkmers);
    REQUIRE(ukm.find(readkmers[0].kmer) == ukm.end()); // FAILS TO FIND NON PRESENT KMERS

    readkmers.clear();
    skf.produce_all_kmers(seqPresent.data(),readkmers);
    REQUIRE(ukm.find(readkmers[0].kmer) != ukm.end()); // FINDS PRESENT KMERS
}

TEST_CASE("UniqueKmerIndex63 create and lookup") {
    std::vector<KmerIDX128> readkmers;
    StreamKmerIDXFactory128 skf(63);
    std::string seqMissing = "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA";
    std::string seqPresent = "CTTGCGGGTTTCCAGGAACTGGCTGTCCTCGGCGTTCAGCGCCATCGACTTCCAGTCCAGCCC";
    WorkSpace ws;
    SequenceDistanceGraph sg(ws);
    sg.load_from_gfa("../tests/datasets/graph/tgraph.gfa");
    Unique63merIndex ukm(sg);

    REQUIRE(ukm.getMap().size() != 0);

    skf.produce_all_kmers(seqMissing.data(),readkmers);
    REQUIRE(ukm.find(readkmers[0].kmer) == ukm.end()); // FAILS TO FIND NON PRESENT KMERS

    readkmers.clear();
    skf.produce_all_kmers(seqPresent.data(),readkmers);
    REQUIRE(ukm.find(readkmers[0].kmer) != ukm.end()); // FINDS PRESENT KMERS
}

TEST_CASE("NKmerIndex create and lookup") {
    const uint8_t k = 15;

    std::vector<uint64_t> readkmers;
    StreamKmerFactory skf(k);
    std::string seqMissing = "AAAAAAAAAAAAAAA";
    std::string seqPresent = "CTTGCGGGTTTCCAG";
    WorkSpace ws;
    SequenceDistanceGraph sg(ws);
    sg.load_from_gfa("../tests/datasets/graph/tgraph.gfa");
    NKmerIndex assembly_kmers(sg,k);

    REQUIRE(!assembly_kmers.empty());

    skf.produce_all_kmers(seqMissing.data(),readkmers);
    REQUIRE(assembly_kmers.find(readkmers[0]) == assembly_kmers.end()); // FAILS TO FIND NON PRESENT KMERS

    readkmers.clear();
    skf.produce_all_kmers(seqPresent.data(),readkmers);
    REQUIRE(assembly_kmers.find(readkmers[0]) != assembly_kmers.end()); // FINDS PRESENT KMERS
}

TEST_CASE("SatKmerIndex create and lookup") {
    const uint8_t k = 10;

    std::vector<uint64_t> readkmers;
    StreamKmerFactory skf(k);
    std::string seqMissing = "AAAAAAAAAA";
    std::string seqPresent = "CTTGCGGGTT";
    WorkSpace ws;
    SequenceDistanceGraph sg(ws);
    sg.load_from_gfa("../tests/datasets/graph/tgraph.gfa");
    SatKmerIndex assembly_kmers(sg,k);

    REQUIRE(!assembly_kmers.contig_offsets.empty());

    skf.produce_all_kmers(seqMissing.data(),readkmers);
    REQUIRE(assembly_kmers.beginCO(readkmers[0]) == assembly_kmers.endCO(readkmers[0])); // FAILS TO FIND NON PRESENT KMERS

    readkmers.clear();
    skf.produce_all_kmers(seqPresent.data(),readkmers);
    REQUIRE(assembly_kmers.beginCO(readkmers[0]) < assembly_kmers.endCO(readkmers[0])); // FINDS PRESENT KMERS
}