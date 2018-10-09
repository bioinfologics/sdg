//
// Created by Luis Yanes (EI) on 09/10/2018.
//

#include <catch.hpp>
#include <sglib/indexers/UniqueKmerIndex.hpp>

TEST_CASE("UniqueKmerIndex create and lookup") {
    unsigned int K(15);
    std::vector<KmerIDX> readkmers;
    StreamKmerIDXFactory skf(K);
    std::string seqMissing = "AAAAAAAAAAAAAAA";
    std::string seqPresent = "CTTGCGGGTTTCCAG";
    SequenceGraph sg;
    sg.load_from_gfa("tests/datasets/tgraph.gfa");
    UniqueKmerIndex ukm(sg, K);

    REQUIRE(ukm.getMap().size() != 0);

    skf.produce_all_kmers(seqMissing.data(),readkmers);
    REQUIRE(ukm.find(readkmers[0].kmer) == ukm.end()); // FAILS TO FIND NON PRESENT KMERS

    readkmers.clear();
    skf.produce_all_kmers(seqPresent.data(),readkmers);
    REQUIRE(ukm.find(readkmers[0].kmer) != ukm.end()); // FINDS PRESENT KMERS
}