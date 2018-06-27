//
// Created by Ben Ward (EI) on 26/01/2018.
//

#include <tuple>
#include "cxxopts.hpp"
#include <sglib/mappers/threader/NodeMapper.h>
#include <sglib/mappers/threader/MappingThreader.h>

int main(int argc, char **argv) {

    // Input options - user specifies the graph file, the reference FASTA and the output prefix.
    std::string graph_filename;
    std::string reference_filename;
    std::string output_dir;
    std::string index_file;
    std::string mappings_file;
    int log_level;
    bool dump_index, dump_unmapped, dump_mappings, dump_paths;
    double unique_kmer_threshold;

    // DEFAULT PARAMETERS... FIXED FOR NOW.

    // TODO: Make these variable.
    uint8_t m(1), n(1), k(31), min_matches(2);

    std::string outdefault("THREADER_1");

// @formatter:off
    cxxopts::Options options("seq-threader", "Locating query sequences in genome assembly graphs.");

    options.add_options()
            ("help", "Print help")
            ("k", "Kmer size", cxxopts::value<uint8_t>(k) -> default_value("31"), "1-31")
            ("r,reference", "Reference FASTA of sequences to locate in graph", cxxopts::value<std::string>(reference_filename), "FASTA - Sequence file")
            ("g,graph", "Genome graph in GFA format", cxxopts::value<std::string>(graph_filename), "GFA file")
            ("o,output", "Output directory name", cxxopts::value<std::string>(output_dir) -> default_value(outdefault))
            ("i", "Dump the unique kmer index file", cxxopts::value<bool>(dump_index))
            ("u", "Dump non-mapped graph nodes", cxxopts::value<bool>(dump_unmapped))
            ("m", "Dump graph node mappings", cxxopts::value<bool>(dump_mappings))
            ("p", "Dump connected paths", cxxopts::value<bool>(dump_paths))
            ("l,logging", "Logging level (0: INFO, 1: WARN, 2: DEBUG)", cxxopts::value<int>(log_level) -> default_value("0"), "log_level")
            ("t,threshold", "% of unique Kmers required to keep mappings", cxxopts::value<double>(unique_kmer_threshold) -> default_value("90"), "unique_kmer_threshold")
            ("mapfile", "File of precomputed mappings", cxxopts::value<std::string>(mappings_file) -> default_value(""))
            ("indexfile", "File containing the unique kmer index for the graph", cxxopts::value<std::string>(index_file) -> default_value(""));
// @formatter:on

    auto result(options.parse(argc, argv));

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

    if (!sglib::check_or_create_directory(output_dir)) {
        std::cout << "Could not find or create output directory: " << output_dir << '.' << std::endl;
        fail = true;
    }

    if (fail) {
        exit(1);
    }

// LOAD GENOME GRAPH...
    SequenceGraph sg;
    sg.load_from_gfa(graph_filename);

// LOAD OR CONSTRUCT A UNIQUE KMER INDEX FOR GRAPH...
    uniqueKmerIndex uki(k);
    if (index_file.empty()) {
        uki.generate_index(sg, k);
        if (dump_index) {
            uki.write_to_file(output_dir + "/unique_kmer_index");
            std::ofstream kmers(output_dir + "/index_kmers.txt");
            uki.write_kmers_to_file(kmers);
            kmers.close();
        }
    } else {
        uki.read_from_file(index_file);
    }

// CONSTRUCT NODE_MAPPER...

    NodeMapper mppr(sg, uki);

// CONSTRUCT SEQUENCE_MAPPER...

// MAP UNIQUE KMERS FROM FASTA SEQUENCES INTO GRAPH, OR READ FROM FILE...
    if (mappings_file.empty()) {
        mppr.mapSequences(reference_filename);
        if (dump_mappings) {
            mppr.write_mappings_to_binary_file(output_dir + "/raw_node_mappings");
        }
    } else {
        mppr.read_mappings_from_binary_file(mappings_file);
    }

    std::ofstream rawprint(output_dir + "/raw_mappings.txt");
    mppr.printMappings(rawprint, true);
    rawprint.close();

    if (dump_unmapped) {
        sglib::OutputLog() << "Dumping unmapped nodes" << std::endl;
        std::ofstream unmapped(output_dir + "/unmapped_nodes_before_filter.txt");
        mppr.print_unmapped_nodes(unmapped);
        unmapped.close();
    }

// FILTER MAPPINGS FOR QUALITY...
    mppr.filterMappings(unique_kmer_threshold);

    std::ofstream filteredprint(output_dir + "/filtered_mappings.txt");
    mppr.printMappings(filteredprint, true);
    filteredprint.close();

    if (dump_unmapped) {
        sglib::OutputLog() << "Dumping unmapped nodes" << std::endl;
        std::ofstream unmapped(output_dir + "/unmapped_nodes_after_filter.txt");
        mppr.print_unmapped_nodes(unmapped);
        unmapped.close();
    }

// CONNECT MAPPINGS INTO PATHS...
    MappingThreader tdr { mppr };

    tdr.thread_mappings();

    std::ofstream graph_paths_out(output_dir + "/mapped_paths_graph.fasta");
    tdr.graph_threads_to_fasta(graph_paths_out, true);
    graph_paths_out.close();

    std::ofstream query_paths_out(output_dir + "/mapped_paths_query.fasta");
    tdr.query_threads_to_fasta(query_paths_out, true);
    query_paths_out.close();

// Calculate reference inclusion by query...
    sglib::OutputLog() << "Finished threading " << reference_filename << " through graph " << graph_filename << '.' << std::endl;
    tdr.calculate_reference_inclusion();

    std::ofstream thread_lengths_out(output_dir + "/thread_lengths.txt");
    tdr.dump_thread_lengths(thread_lengths_out);
    thread_lengths_out.close();
}