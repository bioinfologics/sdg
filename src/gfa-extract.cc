//
// Created by yanesl on 2/16/18.
//



#include <string>
#include <sglib/SequenceGraph.h>
#include <cxxopts.hpp>
#include <sglib/filesystem/check_or_create_directory.h>
#include <sglib/logger/OutputLog.h>

int main(int argc, char * argv[]) {
    std::string gfa_filename, output_prefix;
    std::string dump_mapped, load_mapped;
    unsigned int log_level(4), size_limit(1000), edge_limit(10);
    std::vector<sgNodeID_t> nodes;
    uint64_t max_mem_gb(4);
    bool stats_only(false);
    try {
//@formatter:off
        cxxopts::Options options("gfa-extract", "Extract region from graph");
        options.add_options()
                ("help", "Print help", cxxopts::value<std::string>(), "")
                ("g,gfa", "input gfa file", cxxopts::value<std::string>(gfa_filename), "filepath")
                ("o,output", "output file prefix", cxxopts::value<std::string>(output_prefix), "path")
                ("log_level", "output log level", cxxopts::value<unsigned int>(log_level), "uint")
                ("s,size_limit", "size limit for region to depth_first_search", cxxopts::value<unsigned int>(size_limit), "uint")
                ("e,edge_limit", "limit of edges to depth_first_search", cxxopts::value<unsigned int>(edge_limit), "uint")
                ("n,nodes", "use the following nodes as seeds", cxxopts::value<std::vector<sgNodeID_t>>(nodes), "list of nodes (n1,n2)");
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
        if (result.count("n") != 1) {
            throw cxxopts::OptionException(" please specify the query nodes using the -n --nodes flag");
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
    std::cout << "Starting DFS" << std::endl;
    std::set<sgNodeID_t> subnodes;
    for (const auto &n:nodes) {
        // Go forward on n
        auto v = sg.depth_first_search(n, size_limit, edge_limit, subnodes);
        subnodes.insert(v.cbegin(),v.cend());
        // Go backwards on n
        v = sg.depth_first_search(-n, size_limit, edge_limit, subnodes);
        subnodes.insert(v.cbegin(),v.cend());
    }
    std::copy(subnodes.begin(), subnodes.end(), std::ostream_iterator<sgNodeID_t>(std::cout, ", "));
    std::cout << std::endl;
    std::cout << subnodes.size() << " nodes in solution\n";
    SequenceSubGraph ssg(sg, std::vector<sgNodeID_t>(subnodes.begin(),subnodes.end()));

}