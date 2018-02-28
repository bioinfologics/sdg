//
// Created by Luis Yanes (EI) on 08/12/2017.
//
#include "gtest/gtest.h"
#include <sglib/SMR.h>
#include <sglib/readers/FileReader.h>
#include <sglib/readers/SequenceGraphReader.h>
#include <sglib/factories/KMerIDXFactory.h>
#include <sglib/SequenceGraph.h>

TEST(basic_check, test_smr_from_fasta){
    SMR<KmerIDX,
    kmerIDXFactory<FastaRecord>,
            FastaReader<FastaRecord>,
            FastaRecord, FastxReaderParams, KMerIDXFactoryParams> ref_kmerIDX_SMR({1}, {31}, {4*GB, 0,
                                                                                             100,
                                                                                             "fasta"});
    std::vector<KmerIDX> vector_unique_kmers(ref_kmerIDX_SMR.read_from_file("../tests/datasets/test.fasta"));
    EXPECT_NE(vector_unique_kmers.size(), 0);
}

TEST(basic_check, test_smr_from_fastq){
    SMR<KmerIDX,
            kmerIDXFactory<FastqRecord>,
            FastqReader<FastqRecord>,
            FastqRecord, FastxReaderParams, KMerIDXFactoryParams> ref_kmerIDX_SMR({1}, {31}, {4*GB, 0,
                                                                                             100,
                                                                                             "fastq"});
    std::vector<KmerIDX> vector_unique_kmers(ref_kmerIDX_SMR.read_from_file("../tests/datasets/test.fastq"));
    EXPECT_NE(vector_unique_kmers.size(), 0);
}

TEST(basic_check, test_smr_from_gfa){
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
    EXPECT_NE(unique_kmers.size(), 0);
}
