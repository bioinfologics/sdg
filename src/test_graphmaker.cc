#include <iostream>
#include <fstream>
#include <sglib/WorkSpace.hpp>
#include <sglib/processors/GraphMaker.hpp>
#include "sglib/logger/OutputLog.h"
#include "cxxopts.hpp"


int main(int argc, char * argv[]) {
    std::string seq1="ACGATTTGACACGTAACGATTACGTAACGATTGAGTAACGAACACA";
    std::string seq2="ACGATTTGACACGTAACGATCACGTAACGATTGAGTAACGAACACA";
    std::string seq3="ACGATTTGACACGTAACGATGACGTAACGATTGAGTAACGAACACA";
    std::cout<<"seq1: "<<seq1<<std::endl;
    std::cout<<"seq2: "<<seq2<<std::endl;
    std::cout<<"seq3: "<<seq3<<std::endl;
    SequenceGraph sg1;
    SequenceGraph sg2;
    SequenceGraph sg3;
    uint8_t k=19;
    {
        std::cout << std::endl << "Testing with a single sequence" << std::endl;
        std::vector<uint64_t> kmersv;
        std::unordered_set<uint64_t> kmers;
        StringKMerFactory(seq1, k).create_kmers(kmersv);
        for (auto kmer:kmersv) kmers.insert(kmer);
        SequenceGraph sg;
        GraphMaker gm(sg);
        gm.new_graph_from_kmerset_trivial(kmers, k);
        sg.write_to_gfa("dbg_from_seq1.gfa");
    }
    {
        k=31;
        std::cout << std::endl << "Testing with a single sequence: 128bits" << std::endl;
        std::vector<__uint128_t> kmersv;
        std::unordered_set<__uint128_t> kmers;
        StringKMerFactory128(seq1, k).create_kmers(kmersv);
        for (auto kmer:kmersv) kmers.insert(kmer);
        SequenceGraph sg;
        GraphMaker gm(sg);
        gm.new_graph_from_kmerset_trivial128(kmers,k);
        sg.write_to_gfa("dbg_from_seq1_128.gfa");

    }
    {
        k=33;
        std::cout << std::endl << "Testing with a single sequence: 128bits" << std::endl;
        std::vector<__uint128_t> kmersv;
        std::unordered_set<__uint128_t> kmers;
        StringKMerFactory128(seq1, k).create_kmers(kmersv);
        for (auto kmer:kmersv) kmers.insert(kmer);
        SequenceGraph sg;
        GraphMaker gm(sg);
        gm.new_graph_from_kmerset_trivial128(kmers,k);
        sg.write_to_gfa("dbg_from_seq1_128.gfa");

    }
    {
        std::cout<< std::endl << "Testing with two sequences" << std::endl;
        std::vector<uint64_t> kmersv;
        std::unordered_set<uint64_t> kmers;
        StringKMerFactory(seq1, k).create_kmers(kmersv);
        StringKMerFactory(seq2, k).create_kmers(kmersv);
        for (auto kmer:kmersv) kmers.insert(kmer);
        SequenceGraph sg;
        GraphMaker gm(sg);
        gm.new_graph_from_kmerset_trivial(kmers,k);
        sg.write_to_gfa("dbg_from_seq1_seq2.gfa");
    }
    {
        std::cout<< std::endl << "Testing with three sequences" << std::endl;
        std::vector<uint64_t> kmersv;
        std::unordered_set<uint64_t> kmers;
        StringKMerFactory(seq1, k).create_kmers(kmersv);
        StringKMerFactory(seq2, k).create_kmers(kmersv);
        StringKMerFactory(seq3, k).create_kmers(kmersv);
        for (auto kmer:kmersv) kmers.insert(kmer);
        SequenceGraph sg;
        GraphMaker gm(sg);
        gm.new_graph_from_kmerset_trivial(kmers,k);
        sg.write_to_gfa("dbg_from_seq1_seq2_seq3.gfa");
    }


    return 0;
}

