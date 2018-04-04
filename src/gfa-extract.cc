//
// Created by Luis Yanes on 2/16/18.
//

#include <string>
#include <vector>
#include <sglib/graph/SequenceGraph.hpp>
#include <cxxopts.hpp>
#include <sglib/filesystem/helpers.h>
#include <sglib/logger/OutputLog.h>
#include <sglib/factories/KMerIDXFactory.h>
#include <sglib/factories/ContigBlockFactory.h>

int main(int argc, char * argv[]) {
    std::string help;
    std::string gfa_filename, output_prefix;
    std::string dump_mapped, load_mapped;
    unsigned int size_limit(1000), edge_limit(10);
    unsigned int log_level;
    std::string query_file;
    std::string nodes_input;
    std::string subgraph;
    uint64_t max_mem_gb(4);
    bool stats_only(false);
//@formatter:off
    cxxopts::Options options("gfa-extract", "Extract region from graph");
    options.add_options()
            ("h,help", "Print help", cxxopts::value<std::string>()->implicit_value(""), "")
            ("g,gfa", "input gfa file", cxxopts::value<std::string>(gfa_filename), "file path")
            ("o,output", "output file prefix", cxxopts::value<std::string>(output_prefix), "path")
            ("log_level", "output log level", cxxopts::value<unsigned int>(log_level)->default_value("4"), "uint")
            ("s,size_limit", "size limit in base pairs for region to explore", cxxopts::value<unsigned int>(size_limit)->default_value("1000"), "uint")
            ("e,edge_limit", "limit number of edges to explore", cxxopts::value<unsigned int>(edge_limit)->default_value("10"), "uint")
            ("n,nodes", "use the following node as a seed (this option can be specified multiple times)", cxxopts::value<std::string>(nodes_input), "string")
            ("subgraph", "use the following subgraph as a seed", cxxopts::value<std::string>(subgraph), "file path");
//@formatter:on
    try {
        auto result = options.parse(argc, argv);
        if (0 != result.count("help")) {
            std::cout << options.help({""}) << std::endl;
            exit(0);
        }
        if (result.count("g") != 1) {
            throw cxxopts::OptionException(" please specify graph file using the -g --gfa flag");
        }
        if (result.count("o") != 1) {
            throw cxxopts::OptionException(" please specify output prefix using the -o, --output flag");
        }
        if (result.count("n") == 0 and subgraph.empty()) {
            throw cxxopts::OptionException(" please specify the query nodes using the -n --nodes flag, "
                                                   "a subgraph using the -s --subgraph flag "
                                                   "or a FASTA file using -q --query");
        }
    } catch (const cxxopts::OptionException &e) {
        std::cout << "Error parsing options: " << e.what() << std::endl;
        std::cout << options.help({""}) << std::endl;
        exit(1);
    }

    if (!sglib::check_or_create_directory(output_prefix)) {
        exit(1);
    }

    sglib::OutputLogLevel = sglib::LogLevels::DEBUG;

    sglib::OutputLog() << "Welcome to gfa-extract" << std::endl << std::endl;

    if (gfa_filename.size() <= 4 or gfa_filename.substr(gfa_filename.size() - 4, 4) != ".gfa") {

        throw std::invalid_argument("filename of the gfa input does not end in gfa, it ends in '" +
                                    gfa_filename.substr(gfa_filename.size() - 4, 4) + "'");
    }
    sglib::OutputLogLevel = static_cast<sglib::LogLevels>(log_level);
    SequenceGraph sg;
    sg.load_from_gfa(gfa_filename);
    std::vector<std::string> nodes;

    if (!nodes_input.empty()) {
        std::stringstream nodes_stream(nodes_input);
        std::string node;
        while (std::getline(nodes_stream, node, ',')) {
            nodes.push_back(node);
        }
    }

    if (!subgraph.empty()) {
        SequenceGraph ssg;
        ssg.load_from_gfa(subgraph);
        for (auto n = 1ul; n < ssg.nodes.size(); ++n) {
            nodes.push_back(ssg.oldnames[n]);
        }
    }

    if (!nodes.empty()) {
        auto resultNodes = sg.explore_nodes(nodes, size_limit, edge_limit);

        sglib::OutputLog() << resultNodes.size() << " nodes in solution\n";
        for (const auto &n: resultNodes) {
            sglib::OutputLog(sglib::INFO, false) << sg.oldnames[std::abs(n.node)] << " ";
        }
        sglib::OutputLog(sglib::INFO,false) << std::endl;

        SequenceSubGraph ssg(sg, resultNodes);
        ssg.write_to_gfa(output_prefix+"subgraph.gfa");
    }

    if (!query_file.empty()) {

    }

    sglib::OutputLog() << "Done" << std::endl;
}