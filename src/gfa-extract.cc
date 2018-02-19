//
// Created by yanesl on 2/16/18.
//



#include <string>
#include <sglib/SequenceGraph.h>
#include <cxxopts.hpp>
#include <sglib/filesystem/check_or_create_directory.h>
#include <sglib/logger/OutputLog.h>
#include <sglib/factories/KMerIDXFactory.h>
#include <sglib/readers/SequenceGraphReader.h>
#include <sglib/factories/ContigBlockFactory.h>
#include <sglib/readers/FileReader.h>
#include <sglib/SMR.h>

int main(int argc, char * argv[]) {
    std::string gfa_filename, output_prefix;
    std::string dump_mapped, load_mapped;
    unsigned int size_limit(1000), edge_limit(10);
    unsigned int log_level;
    std::string query_file;
    std::vector<std::string> nodes;
    uint64_t max_mem_gb(4);
    bool stats_only(false);
    try {
//@formatter:off
        cxxopts::Options options("gfa-extract", "Extract region from graph");
        options.add_options()
                ("help", "Print help", cxxopts::value<std::string>(), "")
                ("g,gfa", "input gfa file", cxxopts::value<std::string>(gfa_filename), "filepath")
                ("o,output", "output file prefix", cxxopts::value<std::string>(output_prefix), "path")
                ("log_level", "output log level", cxxopts::value<unsigned int>(log_level)->default_value("4"), "uint")
                ("s,size_limit", "size limit for region to depth_first_search", cxxopts::value<unsigned int>(size_limit)->default_value("1000"), "uint")
                ("e,edge_limit", "limit of edges to depth_first_search", cxxopts::value<unsigned int>(edge_limit)->default_value("10"), "uint")
                ("q,query", "query sequences file", cxxopts::value<std::string>(query_file), "fasta")
                ("n,nodes", "use the following nodes as seeds", cxxopts::value<std::vector<std::string>>(nodes), "list of nodes (n1,n2)");
//@formatter:on
        auto result = options.parse(argc, argv);
        if (result.count("help")) {
            std::cout << options.help({""}) << std::endl;
            exit(0);
        }
        if (result.count("g") != 1) {
            throw cxxopts::OptionException(" please specify graph file using the -g --gfa flag");
        }
        if (result.count("o") != 1) {
            throw cxxopts::OptionException(" please specify output prefix using the -o, --output flag");
        }
        if (result.count("n") == 0 and result.count("q") == 0) {
            throw cxxopts::OptionException(" please specify the query nodes using the -n --nodes flag or a fasta file using -q --query");
        }
        if (result.count("n") > 0 and result.count("q") > 0 ){
            throw cxxopts::OptionException(" the -q,--query flag and -n,--nodes flag are incompatible, please use one or the other");
        }
    } catch (const cxxopts::OptionException &e) {
        std::cout << "Error parsing options: " << e.what() << std::endl << std::endl
                  << "Use option --help to check command line arguments." << std::endl;
        exit(1);
    }
    if (!sglib::check_or_create_directory(output_prefix)) {
        exit(1);
    }

    sglib::OutputLog() << "Welcome to gfa-extract" << std::endl << std::endl;

    sglib::OutputLog() << "Welcome to map-lr" << std::endl << std::endl;

    if (gfa_filename.size() <= 4 or gfa_filename.substr(gfa_filename.size() - 4, 4) != ".gfa") {

        throw std::invalid_argument("filename of the gfa input does not end in gfa, it ends in '" +
                                    gfa_filename.substr(gfa_filename.size() - 4, 4) + "'");
    }
    sglib::OutputLogLevel = static_cast<sglib::LogLevels>(log_level);
    SequenceGraph sg;
    sg.load_from_gfa(gfa_filename);

    std::cout << std::endl;
    std::cout << std::endl;
    std::cout << std::endl;
    if (!nodes.empty()) {
        std::cout << "Starting DFS" << std::endl;
        std::set<sgNodeID_t> subnodes;
        for (const auto &n:nodes) {
            auto id = sg.oldnames_to_ids[n];
            sglib::OutputLog() << "Processing " << sg.oldnames_to_ids[n] << std::endl;
            // Go forward on n
            auto v = sg.depth_first_search(id, size_limit, edge_limit, subnodes);
            subnodes.insert(v.cbegin(), v.cend());
            // Go backwards on n
            v = sg.depth_first_search(-id, size_limit, edge_limit, subnodes);
            subnodes.insert(v.cbegin(), v.cend());
            sglib::OutputLog() << "Done" << std::endl;
        }
        for (const auto &n:subnodes) {
            std::cout << sg.oldnames[std::abs(n)] << " ";
        }
        std::cout << std::endl;
        std::cout << subnodes.size() << " nodes in solution\n";
        SequenceSubGraph ssg(sg, std::vector<sgNodeID_t>(subnodes.begin(), subnodes.end()));
    }
    if (!query_file.empty()) {
        std::cout << "Mapping sequences " << std::endl;
        uint8_t k(13);
        auto max_coverage(1);
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
        uint32_t min_kmers_to_call_match(2);
        uint32_t min_seen_contig_to_write_output(1);
        std::unordered_set<KmerIDX> us(unique_kmers.begin(), unique_kmers.end());
        SMR<Block,
                ContigBlockFactory<FastaRecord>,
                FastaReader<FastaRecord>,
                FastaRecord, FastxReaderParams, FilterSetParams> contigBlock(
                {0},   // ReaderParams
                {output_prefix, k, us, min_kmers_to_call_match, min_seen_contig_to_write_output},
                {maxmem, 1, 100000, output_prefix});

        std::vector<Block> blocks = contigBlock.read_from_file(query_file);
        std::cout << "Num matched block = " << blocks.size();
        std::copy(blocks.begin(),blocks.end(), std::ostream_iterator<Block>(std::cout, " "));
    }
}