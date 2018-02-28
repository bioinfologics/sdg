//
// Created by Luis Yanes (EI) on 26/01/2018.
//

#include <iostream>
#include <fstream>
#include <sglib/factories/ContigBlockFactory.h>
#include <sglib/filesystem/check_or_create_directory.h>
#include <sglib/readers/FileReader.h>
#include <sglib/readers/SequenceGraphReader.h>
#include <sglib/SMR.h>
#include "cxxopts.hpp"


int main(int argc, char * argv[]) {
    std::string gfa_filename, bubble_contigs_filename, output_prefix, long_reads;
    std::string dump_mapped, load_mapped;
    unsigned int log_level(4);
    uint64_t max_mem_gb(4);
    bool stats_only(false);
//@formatter:off
    cxxopts::Options options("map-lr", "LongRead Mapper");
    options.add_options()
            ("help", "Print help", cxxopts::value<std::string>(),"")
            ("g,gfa", "input gfa file", cxxopts::value<std::string>(gfa_filename), "filepath")
            ("o,output", "output file prefix", cxxopts::value<std::string>(output_prefix), "path")
            ("log_level", "output log level", cxxopts::value<unsigned int>(log_level), "uint");
    options.add_options("Long Read Options")
            ("r,long_reads", "input long reads", cxxopts::value<std::string>(long_reads), "filepath")
            ("d,dump_to","dump mapped reads to file",cxxopts::value<std::string>(dump_mapped), "filepath")
            ("l,load_from", "load mapped reads from file", cxxopts::value<std::string>(load_mapped), "filepath")
            ("max_mem", "maximum_memory when mapping (GB, default: 4)", cxxopts::value<uint64_t>(max_mem_gb)->default_value("4"), "GB");
//@formatter:on
    try {
        auto result = options.parse(argc, argv);

        if (result.count("help")) {
            std::cout << options.help({""}) << std::endl;
            exit(0);
        }

        if (result.count("g") != 1 or result.count("o") != 1 ) {
            throw cxxopts::OptionException(" please specify input files and output prefix");
        }
        if (long_reads.empty()) {
            throw cxxopts::OptionException(" please specify a long reads file");

        }
    } catch (const cxxopts::OptionException &e) {
        std::cout << "Error parsing options: " << e.what() << std::endl;
        std::cout << options.help({""}) << std::endl;
        exit(1);
    }
    if (!sglib::check_or_create_directory(output_prefix)) {
        exit(1);
    }

    sglib::OutputLog() << "Welcome to map-lr" << std::endl << std::endl;
    if (gfa_filename.size() <= 4 or gfa_filename.substr(gfa_filename.size() - 4, 4) != ".gfa") {

        throw std::invalid_argument("filename of the gfa input does not end in gfa, it ends in '" +
                                    gfa_filename.substr(gfa_filename.size() - 4, 4) + "'");
    }
    sglib::OutputLogLevel = static_cast<sglib::LogLevels>(log_level);
    max_mem_gb *= GB;
    SequenceGraph sg;
    sg.load_from_gfa(gfa_filename);

    uint8_t k = 15;
    uint16_t min_matches = 10;

    std::cout << "Mapping sequences " << std::endl;
    unsigned int max_coverage(1);
    auto maxmem(4*GB);
    /*
     * Reads the reference/contigs file and generates the set of kmers and which entry number on the fasta they belong to
     * along with their orientation in the form of a int32_t with + for fwd and - for rev, the output is filtered
     * by min < coverage <= max
     */
    SMR<KmerIDX,
            kmerIDXFactory<FastaRecord>,
            GraphNodeReader<FastaRecord>,
            FastaRecord, GraphNodeReaderParams, KMerIDXFactoryParams> kmerIDX_SMR({1,sg}, {k}, {maxmem, 0, max_coverage,
                                                                                                output_prefix});

    std::vector<KmerIDX> unique_kmers;

    // Get the unique_kmers from the file
    unique_kmers = kmerIDX_SMR.process_from_memory();

    std::vector<uint64_t > uniqKmer_statistics(kmerIDX_SMR.summaryStatistics());
    std::cout << "Number of " << int(k) << "-kmers seen in assembly " << uniqKmer_statistics[0] << std::endl;
    std::cout << "Number of " << int(k) << "-kmers that appear more than 0 < kmer_coverage <= " << max_coverage
              << " in assembly " << uniqKmer_statistics[1] << std::endl;
    std::cout << "Number of contigs from the assembly " << uniqKmer_statistics[2] << std::endl;

    /*
     * Instantiate a ContigLinkFactory that takes the unique_kmers as a parameter, a fastq file with long/linked reads
     * and generates a list of links between contigs.
     */

    uint32_t min_kmers_to_call_match(5);
    uint32_t min_seen_contig_to_write_output(1);
    std::unordered_set<KmerIDX> us(unique_kmers.begin(), unique_kmers.end());

    SMR<Block,
            ContigBlockFactory<FastqRecord>,
            FastqReader<FastqRecord>,
            FastqRecord, FastxReaderParams, FilterSetParams> contigBlock(
            {0},   // ReaderParams
            {output_prefix, k, us, min_kmers_to_call_match, min_seen_contig_to_write_output},
            {maxmem, 1, 100000, output_prefix});

    std::vector<Block> blocks = contigBlock.read_from_file(long_reads, true);
    auto blockReaderStats(contigBlock.summaryStatistics());
    std::cout << "Total records generated " << blockReaderStats[0] << std::endl;
    std::cout << "Total filtered records " << blockReaderStats[1] << std::endl;
    std::cout << "Total read sequences " << blockReaderStats[2] << std::endl;
    std::cout << "Total reader filtered sequences " << blockReaderStats[3] << std::endl;

}
