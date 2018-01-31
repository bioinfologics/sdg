//
// Created by Ben Ward (EI) on 26/01/2018.
//

#include "cxxopts.hpp"
#include "sglib/SequenceMapper.h"
#include "sglib/filesystem/check_or_create_directory.h"

int main(int argc, char **argv) {

    // Input options - user specifies the graph file, the reference FASTA and the output prefix.
    std::string graph_filename;
    std::string reference_filename;
    std::string output_prefix;

    // DEFAULT PARAMETERS... FIXED FOR NOW.

    // TODO: Make these variable.
    uint8_t m(1), n(1), k(31), min_matches(2);

    std::string outdefault("REF_THREADER_1");

    cxxopts::Options options("ref-aligner", "Locating reference sequence in genome assembly graphs.");

    options.add_options()("help", "Print help")
    ("r,reference", "Reference FASTA to locate in graph",cxxopts::value<std::string>(reference_filename), "FASTA - Sequence file")
            ("g,graph", "Genome graph in GFA format",cxxopts::value<std::string>(graph_filename), "GFA file")
            ("o,output", "Output file prefix", cxxopts::value<std::string>(output_prefix)->default_value(outdefault),"prefix_dir");

    options.add_options("Skip-mer shape (m every n, total k)")
            ("m,used_bases", "m (1)", cxxopts::value<uint8_t>(m))
            ("n,skipped_bases", "n (1)", cxxopts::value<uint8_t>(n))
            ("k,total_bases", "k (31)",cxxopts::value<uint8_t>(k));

    auto result (options.parse(argc, argv));

    if (0 != result.count("help")) {
        std::cout << options.help({""}) << std::endl;
        exit(0);
    }

    auto fail = false;

    if (graph_filename.empty()) {
        std::cout << "Error: The assembly graph file parameter wasn't specified, " << std::endl
                  << "Use option --help to check command line arguments." << std::endl;
        fail = true;
    } else if (!std::ifstream(reference_filename)) {
        std::cout << reference_filename << " does not exist." << std::endl;
        fail = true;
    }

    if (fail) {
        exit(1);
    }

    if (!sglib::check_or_create_directory(output_prefix)) {
        std::cout << "Could not find or create output directory: " << output_prefix << '.' << std::endl;
        exit(1);
    }


// LOAD GENOME GRAPH...
    SequenceGraph sg;
    sg.load_from_gfa(graph_filename);


// CONSTRUCT SEQUENCE_MAPPER...

    SequenceMapper mppr(sg);
    mppr.map_sequences(1, reference_filename);

    //mppr.print_mappings();

}