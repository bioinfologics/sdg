//
// Created by Luis Yanes on 2/16/18.
//

#include <string>
#include <sglib/SequenceGraph.h>
#include <cxxopts.hpp>
#include <sglib/filesystem/check_or_create_directory.h>
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
    std::vector<std::string> nodes;
    uint64_t max_mem_gb(4);
    bool stats_only(false);
//@formatter:off
    cxxopts::Options options("gfa-extract", "Extract region from graph");
    options.add_options()
            ("h,help", "Print help", cxxopts::value<std::string>()->implicit_value(""), "")
            ("g,gfa", "input gfa file", cxxopts::value<std::string>(gfa_filename), "file path")
            ("o,output", "output file prefix", cxxopts::value<std::string>(output_prefix), "path")
            ("log_level", "output log level", cxxopts::value<unsigned int>(log_level)->default_value("4"), "uint")
            ("s,size_limit", "size limit for region to depth_first_search", cxxopts::value<unsigned int>(size_limit)->default_value("1000"), "uint")
            ("e,edge_limit", "limit of edges to depth_first_search", cxxopts::value<unsigned int>(edge_limit)->default_value("10"), "uint")
            ("n,nodes", "use the following nodes as seeds", cxxopts::value<std::vector<std::string>>(nodes), "list of nodes (n1,n2)");
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
        if (result.count("n") == 0) {
            throw cxxopts::OptionException(" please specify the query nodes using the -n --nodes flag or a FASTA file using -q --query");
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
        std::set<sgNodeID_t> forward_subnodes;
        std::set<sgNodeID_t> backward_subnodes;
        for (const auto &n:nodes) {
            auto id = sg.oldnames_to_ids[n];
            sglib::OutputLog() << "Processing " << sg.oldnames_to_ids[n] << std::endl;
            // Go forward on n
            auto v = sg.depth_first_search(id, size_limit, edge_limit, forward_subnodes);
            forward_subnodes.insert(v.cbegin(), v.cend());
            // Go backwards on n
            v = sg.depth_first_search(-id, size_limit, edge_limit, backward_subnodes);
            backward_subnodes.insert(v.cbegin(), v.cend());
            sglib::OutputLog() << "Done" << std::endl;
        }
        for (const auto &n:forward_subnodes) {
            std::cout << sg.oldnames[std::abs(n)] << " ";
        }
        for (const auto &n:backward_subnodes) {
            std::cout << sg.oldnames[std::abs(n)] << " ";
        }
        std::cout << std::endl;
        std::set<sgNodeID_t> subnodes;
        subnodes.insert(backward_subnodes.begin(),backward_subnodes.end());
        subnodes.insert(forward_subnodes.begin(),forward_subnodes.end());
        std::cout << backward_subnodes.size() << " nodes in solution\n";
        SequenceSubGraph ssg(sg, std::vector<sgNodeID_t>(subnodes.begin(), subnodes.end()));
    }
}