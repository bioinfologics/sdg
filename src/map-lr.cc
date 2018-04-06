//
// Created by Luis Yanes (EI) on 26/01/2018.
//

#include <iostream>
#include <fstream>
#include <sglib/factories/ContigBlockFactory.h>
#include <sglib/filesystem/helpers.h>
#include <sglib/readers/FileReader.h>
#include <sglib/readers/SequenceGraphReader.h>
#include <sglib/SMR.h>
#include <sglib/mappers/LongReadMapper.hpp>
#include "cxxopts.hpp"

#include <sglib/utilities/omp_safe.hpp>
#include <vector>
#include <sglib/mappers/minimap2/minimap.h>

void map_reads(mm_mapopt_t *opt, mm_idx_t *mi, LongReadsDatastore &datastore, std::unordered_set<uint32_t> readIDs) {
    {
        std::ofstream matchHeader("header.paf");
        matchHeader << "read_id\t"
                    << "read_name\t"
                    << "read_len\t"
                    << "read_nmap\t"
                    << "node_id\t"
                    << "node_name\t"
                    << "node_len\t"
                    << "node_start\t"
                    << "node_end\t"
                    << "read_start\t"
                    << "read_end\t"
                    << "strand\t"
                    << "nmatches\t"
                    << "blocklen\t"
                    << "mapq\n";
    }
#pragma omp parallel
{
    mm_tbuf_t *buf = mm_tbuf_init();
    std::ofstream matchOutput(std::string("thread_")+std::to_string(omp_get_thread_num())+std::string(".paf"));
#pragma omp for
        for (uint32_t readID = 1; readID < datastore.size(); ++readID) {
            std::string read_seq(datastore.get_read_sequence(readID));
            std::string read_name(std::to_string(readID));
            int read_len(static_cast<int>(read_seq.length()));
            int n_regs0;
            mm_reg1_t *regs0 = mm_map(mi, read_len, read_seq.data(), &n_regs0, buf, opt, read_name.data());
            for (int j = 0; j < n_regs0; ++j) {
                matchOutput << readID
                        << "\t" << read_name
                        << "\t" << read_len
                        << "\t" << j + 1
                        << "\t" << regs0[j].rid
                        << "\t" << mi->seq[regs0[j].rid].name
                        << "\t" << mi->seq[regs0[j].rid].len
                        << "\t" << regs0[j].rs
                        << "\t" << regs0[j].re
                        << "\t" << regs0[j].qs
                        << "\t" << regs0[j].qe
                        << "\t" << "+-"[regs0[j].rev]
                        << "\t" << regs0[j].mlen
                        << "\t" << regs0[j].blen
                        << "\t" << regs0[j].mapq
                        << "\n";
            }
            for (int i = 0; i<n_regs0;i++) free(regs0[i].p);
            free(regs0);
        }
    mm_tbuf_destroy(buf);
}

}

void map_using_minimap(uint8_t k, SequenceGraph &sg, std::string &long_reads) {
    std::vector<const char *> names(sg.nodes.size());
    std::vector<const char *> seqs(sg.nodes.size());
    for (std::vector<std::string>::size_type i = 1; i < seqs.size(); i++) {
        names[i] = sg.oldnames[i].data();
        seqs[i] = sg.nodes[i].sequence.data();
    }

    mm_mapopt_t opt;
    mm_mapopt_init(&opt);
    mm_idx_t *graph_index = mm_idx_str(15, k, 0, 14, static_cast<int>(seqs.size()-1), &seqs[1], &names[1]);
    mm_mapopt_update(&opt, graph_index);

    LongReadsDatastore datastore(long_reads, "reads_index.idx");
    map_reads(&opt, graph_index, datastore, {});
    mm_idx_destroy(graph_index);
}

void map_using_unique_kmers(uint8_t k, SequenceGraph &sg, std::string &output_prefix, std::string &long_reads){
    uint16_t min_matches = 10;

    std::cout << "Mapping sequences " << std::endl;
    unsigned int max_coverage(1);
    auto maxmem(4 * GB);
    /*
     * Reads the reference/contigs file and generates the set of kmers and which entry number on the fasta they belong to
     * along with their orientation in the form of a int32_t with + for fwd and - for rev, the output is filtered
     * by min < coverage <= max
     */
    SMR<KmerIDX,
            kmerIDXFactory<FastaRecord>,
            GraphNodeReader<FastaRecord>,
            FastaRecord, GraphNodeReaderParams, KMerIDXFactoryParams> kmerIDX_SMR({1, sg}, {k},
                                                                                  {maxmem, 0, max_coverage,
                                                                                   output_prefix});

    std::vector<KmerIDX> unique_kmers;

    // Get the unique_kmers from the file
    unique_kmers = kmerIDX_SMR.process_from_memory(false);

    std::vector<uint64_t> uniqKmer_statistics(kmerIDX_SMR.summaryStatistics());
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

    SMR<ContigBlockFactory<FastqRecord>::Block,
            ContigBlockFactory<FastqRecord>,
            FastqReader<FastqRecord>,
            FastqRecord, FastxReaderParams, FilterSetParams> contigBlock(
            {0},   // ReaderParams
            {output_prefix, k, us, min_kmers_to_call_match, min_seen_contig_to_write_output},
            {maxmem, 1, 100000, output_prefix});

    std::vector<ContigBlockFactory<FastqRecord>::Block> blocks = contigBlock.read_from_file(long_reads, true);
    auto blockReaderStats(contigBlock.summaryStatistics());
    std::cout << "Total records generated " << blockReaderStats[0] << std::endl;
    std::cout << "Total filtered records " << blockReaderStats[1] << std::endl;
    std::cout << "Total read sequences " << blockReaderStats[2] << std::endl;
    std::cout << "Total reader filtered sequences " << blockReaderStats[3] << std::endl;

}

void map_using_sketches(uint8_t k, SequenceGraph &sg, std::string &output_prefix, std::string &long_reads) {
    uint8_t w = 1;
    LongReadsDatastore datastore(long_reads, "reads_index.idx");
    LongReadMapper rm(k, w, sg, datastore);
    rm.map_reads();
}

int main(int argc, char * argv[]) {
    sglib::OutputLog(false) << "Welcome to LongRead Mapper"<<std::endl<<std::endl;
    sglib::OutputLog(false) << "Git origin: " << GIT_ORIGIN_URL << " -> "  << GIT_BRANCH << std::endl;
    sglib::OutputLog(false) << "Git commit: " << GIT_COMMIT_HASH << std::endl<<std::endl;
    sglib::OutputLog() << "Executed command:"<<std::endl;
    for (auto i=0;i<argc;i++) sglib::OutputLog(false) <<argv[i]<<" ";
    sglib::OutputLog(false) <<std::endl<<std::endl;

    std::string gfa_filename, bubble_contigs_filename, output_prefix, long_reads;
    std::string dump_mapped, load_mapped;
    unsigned int log_level(0);
    uint64_t max_mem_gb(4);
    bool stats_only(false);
    uint8_t K(31);
//@formatter:off
    cxxopts::Options options("map-lr", "LongRead Mapper");
    options.add_options()
            ("help", "Print help", cxxopts::value<std::string>(),"")
            ("k,mer_size", "K-mer size for indexing/mapping", cxxopts::value<uint8_t>(K)->default_value("31"), "0-31")
            ("g,gfa", "input gfa file", cxxopts::value<std::string>(gfa_filename), "filepath")
            ("o,output", "output file prefix", cxxopts::value<std::string>(output_prefix), "path")
            ("l,log_level", "output log level", cxxopts::value<unsigned int>(log_level)->default_value("0"), "uint")
            ("r,long_reads", "input long reads", cxxopts::value<std::string>(long_reads), "filepath")
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

    if (gfa_filename.size() <= 4 or gfa_filename.substr(gfa_filename.size() - 4, 4) != ".gfa") {

        throw std::invalid_argument("filename of the gfa input does not end in gfa, it ends in '" +
                                    gfa_filename.substr(gfa_filename.size() - 4, 4) + "'");
    }
    sglib::OutputLogLevel = static_cast<sglib::LogLevels>(log_level);
    max_mem_gb *= GB;
    SequenceGraph sg;
    sg.load_from_gfa(gfa_filename);

    if (0) {
        map_using_minimap(K, sg, long_reads);
    }
    if (1) {
        map_using_sketches(K, sg, output_prefix, long_reads);
    }

    sglib::OutputLog() << "Done!" << std::endl;
    return 0;
}
