//
// Created by Luis Yanes (EI) on 23/08/2018.
//

#include <catch.hpp>
#include <sdglib/readers/FileReader.hpp>
#include <sdglib/workspace/WorkSpace.hpp>
#include <random>

TEST_CASE("Workspace create, read, write") {
    // Long reads
    std::string Lr_filepath("../tests/datasets/workspace/long_reads/long_reads.fastq");
    std::string Lrds_output_path("long_reads.loseq");
    LongReadsDatastore Lrds(Lr_filepath, Lrds_output_path);
    // Linked reads
    std::string lr1_filepath("../tests/datasets/workspace/10x/10x_R1.fastq");
    std::string lr2_filepath("../tests/datasets/workspace/10x/10x_R2.fastq");
    std::string lrds_output_path("10x.lseq");
    LinkedReadsDatastore lrds(lr1_filepath, lr2_filepath, lrds_output_path, LinkedReadsFormat::seq);
    // Paired reads
    std::string r1_filepath("../tests/datasets/workspace/pe/pe_R1.fastq");
    std::string r2_filepath("../tests/datasets/workspace/pe/pe_R2.fastq");
    std::string prds_output_path("pe.prseq");
    PairedReadsDatastore prds(r1_filepath, r2_filepath, prds_output_path);

    WorkSpace out, in;
    out.sdg.load_from_gfa("../tests/datasets/tgraph.gfa");
    out.kci.add_counts_from_datastore(prds);
    out.long_read_datastores.emplace_back(Lrds);
    out.long_read_mappers.emplace_back(out, out.long_read_datastores.back());

    out.paired_read_datastores.emplace_back(prds);
    out.paired_read_mappers.emplace_back(out, out.paired_read_datastores.back(),
            out.uniqueKmerIndex, out.unique63merIndex);

    out.linked_read_datastores.emplace_back(lrds);
    out.linked_read_mappers.emplace_back(out, out.linked_read_datastores.back(),
            out.uniqueKmerIndex, out.unique63merIndex);

    out.dump_to_disk("workspace.bsgws");

    in.load_from_disk("workspace.bsgws");


    REQUIRE( out.sdg == in.sdg);
}

TEST_CASE("Long reads datastore create, read, write") {
    {
        std::string lr_filepath("../tests/datasets/workspace/long_reads/long_reads.fastq");
        std::string lrds_output_path("long_reads.loseq");
        LongReadsDatastore ds(lr_filepath, lrds_output_path);
    }

    LongReadsDatastore ds("long_reads.loseq");
    BufferedSequenceGetter bufferedSequenceGetter(ds);

    //random number engine
    std::mt19937 gen(10); // Always using same seed to get same results
    std::uniform_int_distribution<> dis(1, ds.size());

    std::array<uint64_t, 50> reads_to_check;
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
}

TEST_CASE("10x reads datastore create, read, write") {
    {
        std::string r1_filepath("../tests/datasets/workspace/10x/10x_R1.fastq");
        std::string r2_filepath("../tests/datasets/workspace/10x/10x_R2.fastq");
        std::string lrds_output_path("10x.lseq");
        LinkedReadsDatastore ds(r1_filepath, r2_filepath, lrds_output_path, LinkedReadsFormat::seq);
    }

    const LinkedReadsDatastore ds("10x.lseq");
    BufferedLRSequenceGetter bufferedSequenceGetter(ds, 128*1024,260);

    //random number engine
    std::mt19937 gen(10); // Always using same seed to get same results
    std::uniform_int_distribution<> dis(1, ds.size());

    std::array<uint64_t, 50> reads_to_check;
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
}

TEST_CASE("PE reads datastore create, read, write") {
    {
        std::string r1_filepath("../tests/datasets/workspace/pe/pe_R1.fastq");
        std::string r2_filepath("../tests/datasets/workspace/pe/pe_R2.fastq");
        std::string lrds_output_path("pe.prseq");
        PairedReadsDatastore ds(r1_filepath, r2_filepath, lrds_output_path);
    }

    const PairedReadsDatastore ds("pe.prseq");
    BufferedPairedSequenceGetter bufferedSequenceGetter(ds, 128*1024,260);

    //random number engine
    std::mt19937 gen(10); // Always using same seed to get same results
    std::uniform_int_distribution<> dis(1, ds.size());

    std::array<uint64_t, 50> reads_to_check;
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
    SequenceGraph sg;
    sg.load_from_gfa("../tests/datasets/tgraph.gfa");
    REQUIRE(sg.nodes.size() > 1);
}