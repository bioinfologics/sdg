//
// Created by Luis Yanes (EI) on 23/08/2018.
//

#include <catch.hpp>
#include <sglib/factories/KMerIDXFactory.h>

TEST_CASE("StreamKmerIDXFactory generates all kmers") {
    unsigned int K(15);
    std::string sequence = "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA";
    std::vector<KmerIDX> kmers;
    StreamKmerIDXFactory skf(K);

    skf.produce_all_kmers(sequence.data(), kmers);

    REQUIRE(kmers.size() == sequence.size()-K+1);
}