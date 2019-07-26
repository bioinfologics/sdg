//
// Created by Luis Yanes (EI) on 23/08/2018.
//

#include <catch.hpp>
#include <sdglib/readers/FileReader.hpp>
#include <sdglib/workspace/WorkSpace.hpp>
#include <random>
#include <sdglib/datastores/ReadPathsDatastore.hpp>

TEST_CASE("Workspace create, read, write") {
    // Long reads
    std::string Lr_filepath("../tests/datasets/workspace/long_reads/long_reads.fastq");
    std::string Lrds_output_path("long_reads.loseq");
    LongReadsDatastore::build_from_fastq(Lrds_output_path, Lrds_output_path, Lr_filepath);
    // Linked reads
    std::string lr1_filepath("../tests/datasets/workspace/10x/10x_R1.fastq");
    std::string lr2_filepath("../tests/datasets/workspace/10x/10x_R2.fastq");
    std::string lrds_output_path("10x.lseq");
    LinkedReadsDatastore::build_from_fastq(lrds_output_path, lrds_output_path, lr1_filepath, lr2_filepath, LinkedReadsFormat::raw);
    // Paired reads
    std::string r1_filepath("../tests/datasets/workspace/pe/pe_R1.fastq");
    std::string r2_filepath("../tests/datasets/workspace/pe/pe_R2.fastq");
    std::string prds_output_path("pe.prseq");
    PairedReadsDatastore::build_from_fastq(prds_output_path, r1_filepath, r2_filepath, prds_output_path);

    WorkSpace out, in;
    out.sdg.load_from_gfa("../tests/datasets/graph/tgraph.gfa");
    out.long_reads_datastores.emplace_back(out, Lrds_output_path);

    out.paired_reads_datastores.emplace_back(out, prds_output_path);

    out.linked_reads_datastores.emplace_back(out, lrds_output_path);

    out.add_kmer_counter("kctest", 31).add_count("prtest", out.paired_reads_datastores[0]);

    out.add_operation("test","test","test");

    out.dump_to_disk("workspace.sdgws");

    in.load_from_disk("workspace.sdgws");


    REQUIRE( out.kmer_counters == in.kmer_counters);
    REQUIRE( out.journal == in.journal);
    REQUIRE( out.sdg == in.sdg);

    ::unlink("10x.lseq");
    ::unlink("pe.prseq");
    ::unlink("long_reads.loseq");
    ::unlink("workspace.sdgws");
    ::unlink("kctest.count");
}

TEST_CASE("Long reads datastore create, read, write") {
    {
        std::string lr_filepath("../tests/datasets/workspace/long_reads/long_reads.fastq");
        std::string lrds_output_path("long_reads.loseq");
        LongReadsDatastore::build_from_fastq(lrds_output_path, lrds_output_path, lr_filepath);
    }

    WorkSpace ws;
    LongReadsDatastore ds(ws, "long_reads.loseq");
    ReadSequenceBuffer bufferedSequenceGetter(ds);

    std::string read1("AAAAGAGAGGCAAGCCAAAACGGCGCGCCCAAGTGAAAGGCGCCGACGTCCGTCGACTAAGAGCCCCATAGTAGGGGACGCGCCAAGGTACAGGGCGCCGACTAAGACCAACGTGAAGCGCCGACTCAAGCCCACTTTCATCCAATGTGTCCCACGCAAGCCGAGACGCACAACCATGTGTCGATTTCAATGGAATTCAGTATACACTAATGTCCTCATAGCTCTGAAGGGTGCACAGATGAAGATGCAAATAATGCGTTACCTAGGTGATTGCGTGACACCGGACAATTGAAGTAATTCAACCGGTTATATTGCTTATTCCTCTTGCGAAATAACCTTGGGGGACCTCTTTATGAATCTCCACCTGTCAGCAGATGCGACAAATATCTATTGAATGATTGGGAGCGGTTCAATGCCTCCTAATGTTCAGGTAAAATTGGAGCAGACTGATGGTACTATTACAAACTTTTAACATCAACCAACAAGTTGTGCCAATAGATAATGGAGCAGATATTCGGGTACCATCGATCGGATGTCCAGAATTGTCCTTATCTGGATTCCCGGATCAGGCACAGGAGCAACTGTGATACCGGAACAGCAGCAGGATTCCGATCAAGTACCGATGGTCCTGTGCTCCATGGGATGCTTCGGTGAACTTGGCACGCCAGGTTTCAGCCGGTTCCTTACAATCATTGATTACAAATCGGAAGGATCCATGACTAACTCCGGAATGTGGATAAGTCACACGGAGGTACTTGTCCGATGTATTCGCGAGATCTGGCTGTATCTTTGAGTTTTAATTGAAGTTGATTTCCGTACATGGTAAGAAAACAGGGTAAATCCCGTTAATAGAACACAAAAACAACATTGCAAATAATTATCTATATAATCAAATAATGAGCAAGCGAGTGACGAGCATAAGATAAAAATGAAATAATAAGAACAAGAACAGAGAACTGACAAAATATAAAGAATCGGCAACAGTGGTCGTCTTCGGCAAGATTACTAACAGAAAATACAACGCACGCATGCATACATAACATAACATAAAAATATCGGTGAGCAAGGGAGCTCACGCAGGGGTGGGTGTGCATGCATGCATCATCATCAAAAAGTGTGGCTTGTTTTGTGTTTGATTGTTCTTCACAACCACACTCAACTAATTGTAACATACACTGATACTGTCATAATACATACACATATTTCGCTCTACAAATTATAATACTCTGATTCGTATAATATAACTGTAATAAACACATAAAATGAACTTATTCTATTAATTAGATGCTCTGCTAGCTCTACCTACTAACGAACTGCTGACTCACTAAAAGTAGATTTTAGATCAGAACCTGAGATATTCTTTTTGGAATATTGACTTATTAGCACTTAAGCATATTATAAAGCATTCCAGCTTATTTAGGATATGTGTGTGTGTGTATATAACATAAATAGATACAAGGAATGTTTAGATATCCACATCGTGTTTGAATAAGAAAGTTCCAACATATGGAGAAATAGGTTTTAATAGACTGTTATGTATGTCGCTTCCATTACTCCGACTTCTGTATTTTCTTGTCGTCGAGCTTTCTTGTACTCCTCGAGCTGGTTTGATATAGTTAATGAGAATTAAGAAACAGGAGGTCTTAGTGGCAGACGCAGTCTTTTAAAAGCGGGACACACGTGTTTCAAAGAGAATTCTTGATACCGAAAAGCTTAGGTAGCACATTTGTAACGTGGTTCACCCCCTTGATGAACACGGACGTAAAGAAGATATACATAGTCGGCATTCCAATCCCTTGTAGTCACAACTATGTGCGTGACCTGCCGATGCTCGTACGAAAAATGGAAAAGCTGTGGGATATCATCAAGCTTATTGGTCTTCTTCTTCTCCGCCTCCGCTGTCTTCCTCCGTCATCACTCGTCAGATCTCAAAAGAGGGAATGTTCCAGGTATGTTCCATTACGTTCGAAAGATTTACTATTGTGTACTGTTGACGACAGTAGCTGTACTCGACGCACAGACGTGCTGTGTTGTTCCTTTAGCTTTGCATCGAGCTCCGGTTACTAGGCCGGCCTGATAAAATTTATAGTCTAACACACCCACTAGATGGTTGTCTTTAATGAAGCGTGACTGTCCATTCTCTCCGGTAGGGATTTGGAGACTTTCGTGAATGCTTTAGCAGTGTTGCTGCATCATAGGTTCGGCAAATTTTTCATGCATTTTCGGCCTTGATGGTTGATTTACATTTGCTTCTTTTGCCACCAACGAATCCCCACGTTGTGCAGTATTTGACATAGAGCCTCTGGTGGGAAACCGATTGGTCCGAATTCTGTGCAACAAATAGAAAATTCTCGCTATTTCATATGAGTCTCATAAAAGAAACCCATTGCCAGCGAAGACTACTTTTCACAATTTCCTGTAAGCCGTATTGATAACTTTGTGCTCGTCATCTTGAAAATGAAACTTCCAGAAAGTCAATACATGCATAAATATCTCCATTTTCTCCGAGGTTTCTTAAAGTAGGGTCAACTTTTTGGTTTTACGACACTTATGGATCTAGCGTGAGGGAAGCATGATGTCACCTCGGTTCGGAGTGGTAAAAGAGGACTGCTGACAGATGAAAATGCAGTCTTGAGATTTGCGGCATTGAGGCCAATGTTTGGAATATTCCAGAGGTAATATTCTGCATCTTTCGAGGATCCTAAAGTCTCTTTAATTGGTCAAATAGCAGGGTGCGAAAGTCACCTTATCACGACAAATACTTGAACGATGGCTTTCACTCTATTGCACAATCTTTTACTAATTTTCCATTCATACATTTTCTACAATCACGCCGATCCTTTTGAATGTATGCTGCATCTTTGGTCACTTTCGATCTATATCGAAAATTTGAATACCGTGTTCAACATTTATACTGATTGATGCTGGGGCTCCATCTGAGAAAAGAGTAGCAAAAACACAGATGCGTTCGGAGCATCGTGTCCATCACCTCCAGCCTTATCAATGCTGGAATGTAGAATTACTCCGTGAGTTTTGTCTCTTTTTGATAATTGTTCAAATGTTTTGATAAATGAATCATCAATGGTTTGGTCTATCACACTTGCGTCGTTGCTCATGAATTAGGGTACACAATCATCGCTGCCCATTTCATTGTCTCTCGAATCTTCTTCAAACATGGTAAATAACAGCTTGGGCCTGACAGTCTGATGCACCGCAATGTGAGGCGCTATCAACAATAAAACGGCATAAACCTTACCGCTATACTTGCCAGTGCTATTTGTACCTCTTCATCTAGCTCTCTCATCGCCAAGATGTACCTTCATAGAGAGAGGTTTCGTCCGTGATATAGAGTCACACCTCATGTGTTCGTAACGCCACCGTTTGGAACATCAGTCTGGATTGCGATTCCTTTCTGCGTGTCTTTGATTGCTGCCTGTAGTTCTTGTTTTCAATATTTCCTTTCTTATCTATGCATATAAGATGAATGACAGTTTGCGTATATTTGTTATTTATTTTTTGAAAAAGTGAGATGGACAACGAACAGATTCAATTTCAATCTGAGACAATTATGAAGAATAGGAAGAACAACATTGTTTCGTAGAATCCCCACATGCAGGAGAAAAAGGAGACGAGAAAAGACACCTTCCATTTTGTTTTGATGTCGATTCGTTCGAGATCTGCTAGTCGAGTCGACTCGACTCGTAATTGTCGAAGAATGAAATCACCTCGCCGTTGCGTCTGACATGAGAAATTAATGATTATATCCTCGGTCCGAATGAAGTTGAATCTGTTTTGCTCGAAATGTTTTTATATAGGTCTGGTGCTTTCCGAAGAAGCGAATTCTTTCGATCGTCTTTTTTTTTTTCCATCAATATTTAATTCCCGGACCAGCCAGTGTCATCAGTAATCAGGTAGGTTGTCTGTGAATCAATCAATGTTACTGCGATTCATTGCATAGTAAGGACCATTCTCAAGAACTTCCACGTATTTTATAGATGTAGACACGACTTTGGCAACAAGTTTGACTGATCTTTTGATGATGGATGTTCTTACATTTCACTCTATTCCGTACGCCTCTAACTTAGGTGGTGGAACAATAGGGAGACTGGAATATCCTCAAGATTCAGGACAGAACTATACGAAGGCAGAAGCCAGATTCTAGTTGTATGCCCTATATCCCGCGTATTAAGCTCAATAGTAGCTTCAAGGTTACGAAAGAATCACGAAGATTGATAACTGCTAGTATCCGCATGATGGTAGAATAACAAATGAAAATATATGAGAGCAGCGAAATTTTCATTGGTGTCGCTACCAGCAAGTTAGATTGAATGTTTTGCTTGTTTCACATACCGTAAAAGCAGATTATCATATTTGAGTTGCGGGTTGTGGTGTGGCTGTGTGTTGGCTGTAAAGACAGGGAACGTTGTTTTCCGATCAATTTCGACTGGCAAAAATGTTGGTTAATGCATATCGTATATGCCATAATCAATATTCAACCTGTAAGATCGTATCGGTACAGGTTGTGTCGGATTTGGTCTAGGGTCGCGCGCGAAGCGGCTTTAAGTTTACTTCCGGTAAGACACGCATAATTACGAGTGATTCCTGACGAAAAATTTCAAAATCTCATTTTTTTTTTAAATTCTGTCGATTGGGTCAATTGTACGGAATCTTAATTAGGTGTTTAGCGTGCTTGTGTTTGCCGTGTTAATTTATAAGACGATAATAGCTACCGTGCTACTCGGTAGCGCAAATGCGGTGCACGGAGAAAGTACTATACATGGTAGGTACGGTATTGAGTGCGGAAGGATTGTATGCAATAGTCAGGTATCGCAAAGAAAAGTGCAAGAAGATGGGCGCTCCCGCCGCCCTTGGTTTGAAAAAAAAGCAATTGCGCTTGCCTTGGCGTGAGGCAACTAAGGCTGGTCATTGGGAACTGAAAAAAATCTAAATATCTAATACGGTATGTACCGCATCACAAAAAACTTACTGACTTGCCAGCATGTCTGTAATGCTTACGTACATGTAACTCTACAGTACTGGGTACTACTGTAACTGTACTGAACCAGCCAACAGAAAAAAGGAAGGGAGATCGACGGCTATAAACCGGACGAAGCAAGAGCAAGGAAAGATCAATCAAGATTCAAATGACAAATTCAGGGATATATTATTTCCATTATGATGCATATATTAAATAAAATGATTGCATCATACAAGTATATTTATTTTTCAGCTTGACTAATTGAAAAAGAGTATTGTTTAGTCCATAGATCCAATTCTAATGGATTGTTAGTGAAAAACACTGATTTTAGTATAAAATTATAATCATGTATTTGATAATCGATCCAGAATGAATCGTATATGCACCATGTATCGATGGCAATTCTTGTATGTTTAGTTGATCTTTTTTATGGCTATAAAAACGGATGGAATAAGACGTGAGTAAGTAACATGTTATTTTAGACCAGTTAATAGTACCAAGCATTATTAGGTTTGTAAATTGATGCAAAGTTTTAACCATTACATGGATTTACGTTGATACTCCAAAAGATCATTATTGCATATTGTGTTAAGTGCAGAAGTCCTATTTAATATCATTGTTAGTATAGTTAAGTATGATTCTGCTACTGACAGATTGCAACTGCAACAACTATGTAAGTGCATCGCGTCGGCTGCATCACGAATGCGAACAGGAATGATTTACATGCATTATTAAGGGCTTGTACGAGATAGAGAGAAGAAAATGGATTGCTGAATGAAAAAATAGTAATGAAATCTTGAGAGATTAATTACTAGCCTTACATCGATTCGTACGTACAGTAGGTATGGTCGAGATCGAATATAATAATTGCTTGTAAAGATTGAATTTTAATAAATGATTTATTTTTTTGTGGATAAATTCTTTCATATTGTTATGTTATTTAATGTTATGCGTTGAGTATATGTAGTGTGGCTTTTCCTTTGTTTTATTCTGTTTTATTCTTTTTGTGAGTGGCTTTTTTTGCATTTCTTTATGTATTCTATTAATATAGTGTGAAGCGATTCTTTGTTTAATCATTGTGGTCCTTGTTGTTTTTTGGTCGAACTTTACCATCTGATCAATGGTAGTTCTATGTGTGTGTGTGGTGCGTGTGTGATTGATGGTAATTTCATTTGATATTTAGTGAGTGTCGAGGACTTGTTGACTTGAGGCTTTGCCTGTGCTTCCGATCACAACGAGACGGGTGTCATTCTAGTATTCATTCGTAGCGTAATTTTATTTTTCCGCATGTAATGAACGGGCATATGCTGTAGCACAAAACAAAAAATTGTTTAATCTTAATTGATTTAATCCGTTGTTGAATCTGTAAGAGACCACAGACCTAACACGCCGTACCCAATCGTCCCACCATAATCTCAAAAGTCTCAAAGCAGATCAATGAATCATTACACGACGCTCCGCCAAGCCTCAGTATTCGTGCCACCTCGCGTTGTTCCACACCACATCATTTATCAAACAACGACCAACACCACGACGCTCATACTTCGACTCACGCAGCGCGTTGGTCTTCCACACCCATCATAACATTTCTCCTCACCATATCATGTAGTGCATTTCACACCAAGACCATTAAGATAAGGTCTATCTACAACCCTTATTCTTCTTTGTCGGTAGCTCTCAATATATCCATCATTATCATAGCACCACCTTCCTCTATATATATTCCTCATATCGAACTCATTTAGATAATCATTCATGTTTTTGGTTTGATCTTTGAATAACCAACAGCACGAATTTTAGATGAGACACCATTATATACACTGGAGAAATGGAAAATGGGAGAGATAAATTGATAGAGAGAATACCGAGAAACAAATTCGGGTCAGACAAAAACACAGGGAAGATTTTTAATAGACAATGTTTATTTCAATTATTTACCTTTCCTTAGAATTGTTCATCTTGATAACCCACCACACTGACCGAATTCTTCATAGGCATCAGGATCGACGTCCCATGCATATCTGTTAAAGTCACTCCCTACCACGGACCACACCATACGATCACGTGCACCTCACCAAGACGACTTCCATATCCTAAGTCATAGGTCACGGCGGTCTTCCTCTTGTTCAGTTCCAGTCACACGTTCAATTCTTCAGATCACCCATTCTGGAAAAAAATCATTATTATTTGACCACTCCATGCCACCACCCATATACGTCTCCAGATCTCCATATAAGTTGTTGGCGTCCACGTGCACGGATCGACATTTTCACGACTTACTCGATAGAGGCGAATCTTGGCATCGATCGACCATTCTTGGAGGATCATCCAACCATTGATGCCATGCGTGTGGCTCGTGGTGTGTATATTTAACCGGTACCTATTTAAATTACAGCGAATAATCATAATAACTCAGCAATGATATCACACAACGCACTACCTCCATGGACTTCACCGCACACTACCTCTTCGAGCAATATGCTGAATGATGTTGTTGTTTCTAAAACTTAGACTGACAGTGGTGATGATGGATTGTGGTTTTAAAAGCAACAGCGTTGTTAATGGTGATATTATTATGAGTCAGGCAGTAGCAGAGCAATTTAGGCAAGATAGCTGTCGTGTTCTTCAGGACATATATGGCTGTGGATCTTAGTCTTTTTTTCTTTCATGATTTAGAATTTAAATTTAAAGGATAACGCTAGTACTGTATGTCATGAGAGGTTATAATATGTGTGATGTTGTCTTTTTCTTTGCTGTTCGTTCAGATTGAGTTGACTCTCCTCGACTCGTTGCTTCTTCTTGTGTTGGTTGGTGTCTTTCTTTGTTTGTGTGATGCCTTTCTTTTTCTGAGTGGATGAATGGAAATCGAGCAGGTGTTCCTAAACCAGGAAGTGCTTTACGGATATCGAACACGATTTATTACACAAAGCACAAAGGGGAAAGGTACATACAATAATAGGAGGCGCGAAGGTAAACACACGATAATACCGGCACGGAAATACGATACAAAAGTTTTTGGATCTGACTTTTTGATTGGATATAATTTTTGAATATTCAAAAATTGAGTAAGGCTAATAGCAAGTCAGCCAGTAAGAAGAGTCGAGCAGAAGAAGCATAGGGGAAGTGAGTCTTGGAAGTGAGTCGAGTCATGCAGTAATCTTCTGGTTCGTGCTGGCGTGGTAACGCGTCGTATGACTGGTATGTTCAGTACAAGAAGGAATGGTAGAATGATAACTCTAGCGTTTTACGATATGACCATACATCTCTCATATCGAAACTACGACACGATTAACATCTGGTGTGTGTGTGTGTCTGCGTGTATGTGGTTCTTGCCGTGAATCATGTGTTACTATGAAGCAGAAGGAGTACGAAGAGCTTGCGTGTGTGTCAATGAAAAAATATTAAATATTAAAATTGATTTTAACAACACGGAAATACGAAGACAATCTAGATAGTTGTTTTTGTTTTTTCTAAAAAAAGGCATTATATTTTTTTTTAGTAGTCTATTAGTTGCATCTCTGCCGTGAAAGAGATCCAGCATATTTACCTATTTACCAAGTGTATGTGCGCCGGTGTAGTGCATTGTTGCGGTCGTAAGGTGTATGATCCAGTAAAAGCATAGTGGGTGGTGTATGTTATTTACACTACCCATCGATTTCCTTTCTGATTCCATCATAAAGTAAGATGTCTTTTGGTATTATTTATTACACCACATACGCCATCCCTTTTTTTGGAAAGACGACGGATCGAGCTAGCAGCCAACAGTAGTATATCAACATCAGGTCCCTTCTTTTGAGAGTGGTAGTTGACGCAGGGGATGGCACATGCCACAGTATATAAGCGTAATCGCATCGTTGATCTACTCTTTTGTTGTGCTAGATTTAAGGATTCCGGTAATCATCACGTATATCATCAGTACAGATAGTATGTAACCATTATACTACAAATAGTTGAATGTTGTTGAGCGAACTGTCTTAGAGGGTGTTGGTGGTTCGAGATAATAACTTCGTGATGTTTTATATTGGGAATATTCGACGAGTTTAATAGACCAATGATGTTATATTCGGCGGCTGCTTTGCTGCTGCTTGCTTGTTGACTTCCATTCTTTCTTCTTCTGATCATGAATCTTAAAGTTGTTTTTTTGTCAGA");

    std::string ds1(ds.get_read_sequence(1));

    REQUIRE(ds1 == read1);

    std::string read300("AAACCTCGTTAAGAAAATAGTTTTGGGATTGTGTGGTTTTTGGCTTTGGGCCCTTTATGTCTAATTAATTGTATCAAAAAAATAGAACGAAAAATCATATATTTTTCTGACCCTGCCGCGCCCCAAACAATTTTTTTTAAACAAAATCAGTGACAATTCGTGGGGAAGTTGATCTTATGTCATTTAAACATGCAATATTCCATTGAGAAAACCCAACTAACTGGTAAGTTGTACCCTCATCTCCGTCATCTCATGGTATGTAGTAGGTGACAATTTACCGTACGGAATAGAACAAGTAGCATTCAGTTTACGTAGTGGTGCATGAACACTGTAGTACCTTATGCGGCAGTCATTTGTTTATTAGTCTTATCAACACAAACAAACTTAAGCAACAATACCTATGATGAAAGCAGACAATTACTTAAATAAATTCTTGAAAGAAGCAGTTACAGGATAAGCAAAGCATAAAAATATTTGTCTCAACAGCACTACAGTAGCGAATACCAAGTGGATGGCCCTGAAATAGGTAATACTCTAATTTGAAAGCCTTAGACCAATAGTGAATCTCATTATTGCAAACGATCGCTCAGGGCCAGTCTATCTAATAGCTTACACAGCAGAGTGTCTTGTCACAGAATTGTCCTTATAGCCTCTAGTCAACCAAAACTTGGGTACGAATTCAGTTTTGTAGAATGTGATCCAGACCTTAAATAAAATAATCTGGCACCTCTATGCATGCTGGGTTTGAGTCCGTCATTCTCCATTTTCCCCCTTACATTATAAACTACATGGTTAATTTTTCAATTCATCGGATGATCTACTGTTGGTTTTGTTTGGAATGGGAGTGCATTGGTACAAATGTCAGCAACAAGACTAGTTAATTTGGAGAATTGAGGAATGAGAGAACGTTTGTGCAACAGGATAGGACGAAATAATACAGGTGCAACAAATGAAAAATACACTAACAGACCTTGTCGAAATCATGGGTTCGGTTAACAGATTCCAGTTATGGCTCTCAGCTGAAGAGGCCACTTTTTTTGGGGTCAATTTTTTTTGTGGGGTGAATTTCTCAAAAAATGGACTTACACCCCCCTTTGGAGACCAGCGCAGTCTCCGAAGTCTAGGGCTACGGGTCTTTTTATATACTCAATTATGGAAGTCTCGTTGCTTTGTGATGGGCCGCACCTTAACCAGCGTAATACGGGCAAAATTAAATTCGGCTAAAGCCCATAGTAGGGTAATATTCGCGGTGAGCGTTGGGGATTGAACCAGTGCAATGTGTTAAAAGCAGCTGAAATATCTAAGTGAATTTCTGCGGCGATTTTTTTAATTTTTAAAAAACTCAAATAGTTGTTTAATAAAGACTGAGTAAGTAGTCAATGTAATAACTTCGGTATTGTGGTTGAGTTGGCACTGGCCGCCAGCGCGGCTTGGAGATAGGTCCTCGTTTCGATCACTTCCGATGTAATACATATTACTTTGTATAGTATACTTTTAGTTATATCAGGCTTGCAAAGATTTCATTTTTATCTATTTTTAAAGCAAAATACTGCGTGTTATGTAACACGCGTAAAATCTCCGAATAGAATTCCTTATAAAGGAACCTGTTAATGGTTGTCTTCGCATATCTTGAGGAGGAATGGAATGGGTCTCAGCTTGTCATAGTCAATGCATCACGTTCTTATATAATATTCAATGAATTTTTTGTTATTAATCCTACGGAAGCCAGTGCCTTTTGCGGGAAAATTATCCGAAAAAATCATATAAGGATATGAGTAAGGATGAAACAATACTCATACAACAAAATTAATGCCCAAATGCAATGTGAATGGAAATATCTATTTGTGGGATTTATAATATACGGTGTTAGTAGATAATTGATAATAATGTATTGCACTGGTACTGTTTTGCCTTGTTGCTTTGTGTATGCAAATAGCAGTATGTAGACATGATTAGTCATACATTTGATATCAGATTGCCATTTGATTGGGGATCCGGTTTTTTGAATACACCAAGTCCCTCTCCGTAAATATGAGATGGATTGTTGAATCACGCCATTTCTTTGTCGTATTATGTCAACAATTCAAGCTATAGAATTTGGTGTATTACCCGATAAAAAACGCGTTAAGTTGCAGAGAACTATTGCGGGAAAGATGAAAATTCCGAGGGCAAATAGAAGTAGATGACAACATTATGTAAGATTGACACATGAAGAGATAAACCGAAATCCCTCACCTGTGATAAACACAAAGGAGTTTGGAAAATAAACATCATAAATACTTTTGATATTTAGCAAAATAATTATACGCACATATTCTGCACTGTTTCCTCTTGGTTGCTTACACATTCGCTGGGATTCAATTCCCTGCCAACTCTTCCCAAATTTTTCCTAGGGGGCTCTTTCATTCATGTATAAGTTTTTACAGCTTATTACACATAGGTTGCGGCCCCATACAAGAGCAATGAAAGGACATAATATTGAATAATAGAAAACGACACTGCCTTGAACGGATGAGATGGGACGCTAATGTCCGATTTTTGGGGGGAGTAAACATCCCTGG");

    std::string ds300(ds.get_read_sequence(300));

    REQUIRE(ds300 == read300);

    ::unlink("long_reads.loseq");
}

TEST_CASE("10x reads datastore create, read, write") {
    {
        std::string r1_filepath("../tests/datasets/workspace/10x/10x_R1.fastq");
        std::string r2_filepath("../tests/datasets/workspace/10x/10x_R2.fastq");
        std::string lrds_output_path("10x.lseq");
        LinkedReadsDatastore::build_from_fastq(lrds_output_path, lrds_output_path, r1_filepath, r2_filepath, LinkedReadsFormat::raw);
    }

    WorkSpace ws;
    const LinkedReadsDatastore ds(ws, "10x.lseq");
    ReadSequenceBuffer bufferedSequenceGetter(ds, 128*1024,260);

    //random number engine
    std::mt19937 gen(10); // Always using same seed to get same results
    std::uniform_int_distribution<> dis(1, ds.size());

    std::array<uint64_t, 50> reads_to_check{};
    for (auto &r:reads_to_check){
        r = dis(gen);
    }
    const char* seq_ptr;
    std::unordered_map<uint32_t, size_t> read_sz1;
    StringKMerFactory skf(15);
    std::vector<std::pair<bool, uint64_t>> read_kmers;
    for (const auto &r: reads_to_check){
        read_kmers.clear();
        seq_ptr = bufferedSequenceGetter.get_read_sequence(r);
        skf.create_kmers(seq_ptr, read_kmers);
        read_sz1[r]=read_kmers.size();
    }

    bool equal=true;
    uint32_t it=0;
    int i=0;
    for (i = 0; i < reads_to_check.size() and equal; i++) {
        seq_ptr = bufferedSequenceGetter.get_read_sequence(reads_to_check[i]);
        read_kmers.clear();
        skf.create_kmers(seq_ptr, read_kmers);
        if (read_kmers.size() != read_sz1[reads_to_check[i]]) equal = false;
    }

    if (!equal) {
        std::cout << "Failed during iteration " << it << ", on read " << reads_to_check[i];
    }

    REQUIRE(equal);
    ::unlink("10x.lseq");
}

TEST_CASE("PE reads datastore create, read, write") {
    {
        std::string r1_filepath("../tests/datasets/workspace/pe/pe_R1.fastq");
        std::string r2_filepath("../tests/datasets/workspace/pe/pe_R2.fastq");
        std::string prds_output_path("pe.prseq");
        PairedReadsDatastore::build_from_fastq(prds_output_path, r1_filepath, r2_filepath, prds_output_path);
    }

    WorkSpace ws;
    const PairedReadsDatastore ds(ws,"pe.prseq");
    ReadSequenceBuffer bufferedSequenceGetter(ds, 128*1024,260);

    //random number engine
    std::mt19937 gen(10); // Always using same seed to get same results
    std::uniform_int_distribution<> dis(1, ds.size());

    std::array<uint64_t, 50> reads_to_check{};
    for (auto &r:reads_to_check){
        r = dis(gen);
    }
    const char* seq_ptr;
    std::unordered_map<uint32_t, size_t> read_sz1;
    StringKMerFactory skf(15);
    std::vector<std::pair<bool, uint64_t>> read_kmers;
    for (const auto &r: reads_to_check){
        read_kmers.clear();
        seq_ptr = bufferedSequenceGetter.get_read_sequence(r);
        skf.create_kmers(seq_ptr, read_kmers);
        read_sz1[r]=read_kmers.size();
    }

    bool equal=true;
    uint32_t it=0;
    int i=0;
    for (i = 0; i < reads_to_check.size() and equal; i++) {
        seq_ptr = bufferedSequenceGetter.get_read_sequence(reads_to_check[i]);
        read_kmers.clear();
        skf.create_kmers(seq_ptr, read_kmers);
        if (read_kmers.size() != read_sz1[reads_to_check[i]]) equal = false;
    }

    if (!equal) {
        std::cout << "Failed during iteration " << it << ", on read " << reads_to_check[i];
    }

    REQUIRE(equal);
    ::unlink("pe.prseq");
}

TEST_CASE("Fastq file reader") {
    FastqReader<FastqRecord> fastqReader({0}, "../tests/datasets/test.fastq");
    FastqRecord read;
    bool c;
    c = fastqReader.next_record(read);
    uint num_reads(0);
    while (c) {
        num_reads++;
        c = fastqReader.next_record(read);
    }

    REQUIRE(num_reads==10);
}

TEST_CASE("Fasta file reader") {
    FastaReader<FastaRecord> fastaReader({0}, "../tests/datasets/test.fasta");
    FastaRecord read;
    bool c;
    c = fastaReader.next_record(read);
    uint num_reads(0);
    while (c) {
        num_reads++;
        c = fastaReader.next_record(read);
    }

    REQUIRE(num_reads==3);
}

TEST_CASE("Load GFA") {
    sdglib::OutputLogLevel = sdglib::DEBUG;
    WorkSpace ws;
    SequenceDistanceGraph sg(ws);
    sg.load_from_gfa("../tests/datasets/graph/tgraph.gfa");
    REQUIRE(sg.nodes.size() > 1);
}

TEST_CASE("Load GFA2") {
    sdglib::OutputLogLevel = sdglib::DEBUG;
    WorkSpace ws;
    SequenceDistanceGraph sg(ws);
    sg.load_from_gfa("../tests/datasets/graph/test_gfa2.gfa");
    REQUIRE(sg.nodes.size() > 1);
}

TEST_CASE("Load ReadPaths") {

    ReadPathsDatastore ds("../tests/datasets/test.paths");

    REQUIRE(!ds.read_paths.empty());
}