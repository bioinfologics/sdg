//
// Created by Ben Ward (EI) on 26/01/2018.
//

#include <tuple>
#include "cxxopts.hpp"
#include <sglib/SequenceThreader.h>
#include <sglib/filesystem/check_or_create_directory.h>
#include <sglib/logger/OutputLog.h>
#include <sglib/aligners/submat/SubstitutionMatrix.h>
#include <sglib/aligners/algorithms/SmithWaterman.h>

int main(int argc, char **argv) {

    // Input options - user specifies the graph file, the reference FASTA and the output prefix.
    std::string graph_filename;
    std::string reference_filename;
    std::string output_prefix;
    int log_level;
    bool dump_unmapped, dump_mappings, dump_paths;


    // DEFAULT PARAMETERS... FIXED FOR NOW.

    // TODO: Make these variable.
    uint8_t m(1), n(1), k(31), min_matches(2);

    std::string outdefault("THREADER_1");

// @formatter:off
    cxxopts::Options options("seq-threader", "Locating query sequences in genome assembly graphs.");

    options.add_options()
            ("help", "Print help")
            ("k", "Kmer size", cxxopts::value<uint8_t>(k)->default_value("31"), "1-31")
            ("r,reference", "Reference FASTA of sequences to locate in graph", cxxopts::value<std::string>(reference_filename), "FASTA - Sequence file")
            ("g,graph", "Genome graph in GFA format", cxxopts::value<std::string>(graph_filename), "GFA file")
            ("o,output", "Output file prefix", cxxopts::value<std::string>(output_prefix) -> default_value(outdefault), "prefix_dir")
            ("u", "Dump non-mapped graph nodes", cxxopts::value<bool>(dump_unmapped))
            ("m", "Dump graph node mappings", cxxopts::value<bool>(dump_mappings))
            ("p", "Dump connected paths", cxxopts::value<bool>(dump_paths))
            ("l,logging", "Logging level (0: INFO, 1: WARN, 2: DEBUG)", cxxopts::value<int>(log_level) -> default_value("0"), "log_level");
// @formatter:on

    auto result (options.parse(argc, argv));

    if (0 != result.count("help")) {
        std::cout << options.help({""}) << std::endl;
        exit(0);
    }

    sglib::OutputLogLevel = (sglib::LogLevels) log_level;

    auto fail = false;

    if (graph_filename.empty()) {
        std::cout << "Error: The assembly graph file parameter wasn't specified, " << std::endl
                  << "Use option --help to check command line arguments." << std::endl;
        fail = true;
    } else if (!std::ifstream(reference_filename)) {
        std::cout << reference_filename << " does not exist." << std::endl;
        fail = true;
    }

    if (!sglib::check_or_create_directory(output_prefix)) {
        std::cout << "Could not find or create output directory: " << output_prefix << '.' << std::endl;
        fail = true;
    }

    if (fail) {
        exit(1);
    }

// LOAD GENOME GRAPH...
    SequenceGraph sg;
    sg.load_from_gfa(graph_filename);

// CONSTRUCT SEQUENCE_MAPPER...
    SequenceThreader tdr(sg, k);
    // Thread FASTA sequences into graph.
    tdr.map_sequences(reference_filename, output_prefix);

// CONNECT MAPPINGS...
    tdr.thread_mappings();

// TEST THE CONNECT_THREADS METHOD SO FAR...
    tdr.bridge_threads();

// DUMP OUTPUT...

    // Mappings dump.
    if (dump_mappings) {
        sglib::OutputLog(sglib::LogLevels::INFO) << "Dumping all mappings." << std::endl;
        std::ofstream mapping_dump(output_prefix + "/mappings.txt");
        tdr.print_mappings(mapping_dump, true);
        mapping_dump.close();
    }

    // Unmapped node diagnostics dump.
    if (dump_unmapped) {
        sglib::OutputLog(sglib::LogLevels::INFO) << "Dumping unmapped nodes." << std::endl;
        std::ofstream unmapped(output_prefix + "/unmapped_nodes.txt");
        tdr.print_unmapped_nodes(unmapped);
        unmapped.close();
    }

    // Paths dump.
    if (dump_paths) {
        sglib::OutputLog(sglib::LogLevels::INFO) << "Dumping paths." << std::endl;
        std::ofstream pathdumpout(output_prefix + "/paths_list.txt");
        tdr.print_paths(pathdumpout, true);
        pathdumpout.close();
    }

    // Print out mapped FASTA sequences.

    std::ofstream graph_paths_out(output_prefix + "/mapped_paths_graph.fasta");
    tdr.graph_threads_to_fasta(graph_paths_out, true);
    graph_paths_out.close();

    std::ofstream query_paths_out(output_prefix + "/mapped_paths_query.fasta");
    tdr.query_threads_to_fasta(query_paths_out);
    query_paths_out.close();

    std::ofstream bridged_graph_paths_out(output_prefix + "/bridged_mapped_paths_graph.fasta");
    tdr.bridged_graph_threads_to_fasta(bridged_graph_paths_out);
}