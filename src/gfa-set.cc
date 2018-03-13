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
    std::string operation;
    std::string gfa_filename, output_prefix;
    unsigned int size_limit(1000), edge_limit(10);
    unsigned int log_level;
    std::string subgraphA, subgraphB;
    uint64_t max_mem_gb(4);
    bool stats_only(false);
//@formatter:off
    cxxopts::Options options("gfa-set", "Graph set tool");
    options.add_options()
            ("h,help", "Print help", cxxopts::value<std::string>()->implicit_value(""), "")
            ("g,gfa", "input gfa file", cxxopts::value<std::string>(gfa_filename), "file path")
            ("o,output", "output file prefix", cxxopts::value<std::string>(output_prefix), "path")
            ("log_level", "output log level", cxxopts::value<unsigned int>(log_level)->default_value("4"), "uint")
            ("op", "operation", cxxopts::value<std::string>(operation), "[union,intersection,difference]")
            ("a,subgraphA", "subgraph A", cxxopts::value<std::string>(subgraphA), "file path")
            ("b,subgraphB", "subgraph B", cxxopts::value<std::string>(subgraphB), "file path");
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
        if (result.count("a") != 1 or result.count("b") != 1) {
            throw cxxopts::OptionException(" please specify both a and b graphs");
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

    sglib::OutputLog() << "Welcome to gfa-set" << std::endl << std::endl;

    if (gfa_filename.size() <= 4 or gfa_filename.substr(gfa_filename.size() - 4, 4) != ".gfa") {

        throw std::invalid_argument("filename of the gfa input does not end in gfa, it ends in '" +
                                    gfa_filename.substr(gfa_filename.size() - 4, 4) + "'");
    }
    sglib::OutputLogLevel = static_cast<sglib::LogLevels>(log_level);
    SequenceGraph sg;
    sg.load_from_gfa(gfa_filename);

    SequenceGraph a,b;
    a.load_from_gfa(subgraphA);
    b.load_from_gfa(subgraphB);

    std::unordered_set<sgNodeID_t > nodes_a, nodes_b;
    std::unordered_set<sgNodeID_t > links_a, links_b;
    for (const auto &n:a.nodes)
        nodes_a.insert(n);
    for (const auto &n:b.nodes)
        nodes_b.insert(n);

    for (const auto &l:a.links)
        links_a.insert(l);
    for (const auto &l:b.links)
        links_b.insert(l);

    enum string_code {
        union_,
        intersection_,
        difference_,
    };

    string_code hashit (std::string const& inString) {
        if (inString == "union") return union_;
        if (inString == "intersection") return intersection_;
        if (inString == "difference") return difference_;
        std::logic_error("Function not yet implemented");
    }

    std::unordered_set<sgNodeID_t > result_nodes;
    std::unordered_set<Link> result_links;

    switch (hashit(operation)){
        case intersection_:
            std::set_intersection(nodes_a.begin(),nodes_a.end(),nodes_b.begin(), nodes_b.end(), result_nodes);
            std::set_intersection(links_a.begin(),links_a.end(),links_b.begin(), links_b.end(), result_links);
            break;
        case union_:
            std::set_union(nodes_a.begin(),nodes_a.end(),nodes_b.begin(), nodes_b.end(), result_nodes);
            std::set_union(links_a.begin(),links_a.end(),links_b.begin(), links_b.end(), result_links);
            break;
        case difference_:
            std::set_difference(nodes_a.begin(),nodes_a.end(),nodes_b.begin(), nodes_b.end(), result_nodes);
            std::set_difference(links_a.begin(),links_a.end(),links_b.begin(), links_b.end(), result_links);
            break;
    }

    SequenceSubGraph ssg(sg, std::vector<sgNodeID_t >(result_nodes.begin(), result_nodes.end()));
    ssg.write_to_gfa("result.gfa");
    sglib::OutputLog() << "Done" << std::endl;
}