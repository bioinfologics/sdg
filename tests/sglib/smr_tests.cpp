//
// Created by Luis Yanes (EI) on 11/06/2018.
//

#include <catch.hpp>
#include <sdglib/SMR.hpp>
#include <sdglib/types/KmerTypes.hpp>
#include <sdglib/readers/FileReader.hpp>
#include <sdglib/readers/SequenceGraphReader.hpp>
#include <sdglib/factories/KMerIDXFactory.hpp>
#include <sdglib/graph/SequenceDistanceGraph.hpp>

TEST_CASE("Test SMR from fasta"){
    sdglib::OutputLogLevel = sdglib::DEBUG;
    SMR<KmerIDX,
        kmerIDXFactory<FastaRecord>,
        FastaReader<FastaRecord>,
        FastaRecord, FastxReaderParams, KMerIDXFactoryParams> ref_kmerIDX_SMR({1}, {31}, {4*GB, 0,
                                                                                          100,
                                                                                          "fasta"});
    std::vector<KmerIDX> vector_unique_kmers(ref_kmerIDX_SMR.read_from_file("../tests/datasets/test.fasta"));
    REQUIRE(vector_unique_kmers.size() != 0);
}

TEST_CASE("Test SMR from fastq"){
    sdglib::OutputLogLevel = sdglib::DEBUG;
    SMR<KmerIDX,
            kmerIDXFactory<FastqRecord>,
            FastqReader<FastqRecord>,
            FastqRecord, FastxReaderParams, KMerIDXFactoryParams> ref_kmerIDX_SMR({1}, {31}, {4*GB, 0,
                                                                                              100,
                                                                                              "fastq"});
    std::vector<KmerIDX> vector_unique_kmers(ref_kmerIDX_SMR.read_from_file("../tests/datasets/test.fastq"));
    REQUIRE(vector_unique_kmers.size() != 0);
}

TEST_CASE("Test SMR from memory"){
    sdglib::OutputLogLevel = sdglib::DEBUG;
    SequenceGraph sg;
    sg.load_from_gfa("../tests/datasets/tgraph.gfa");
    SMR<KmerIDX,
            kmerIDXFactory<FastaRecord>,
            GraphNodeReader<FastaRecord>,
            FastaRecord, GraphNodeReaderParams, KMerIDXFactoryParams> kmerIDX_SMR({1,sg}, {31}, {4, 0, 1000,
                                                                                                 "sg"});
    std::vector<KmerIDX> unique_kmers;

    // Get the unique_kmers from the file
    unique_kmers = kmerIDX_SMR.process_from_memory();
    REQUIRE(unique_kmers.size() != 0);
}
