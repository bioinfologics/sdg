//
// Created by Luis Yanes on 2/16/18.
//

#include <string>
#include <sglib/SequenceGraph.h>
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
    std::vector<std::string> nodes;
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
            ("n,nodes", "use the following node as a seed (this option can be specified multiple times)", cxxopts::value<std::vector<std::string>>(nodes), "string")
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

    if (nodes.empty() and !subgraph.empty()) {
        SequenceGraph ssg;
        ssg.load_from_gfa(subgraph);
        for (auto n = 1ul; n < ssg.nodes.size(); ++n) {
            nodes.push_back(ssg.oldnames[n]);
        }
    }

    std::cout << std::endl;
    std::cout << std::endl;
    std::cout << std::endl;
    sglib::OutputLog() << "Starting DFS" << std::endl;
    if (!nodes.empty()) {
        std::set<nodeVisitor> resultNodes;
        // For each node in the list
        for (const auto &n:nodes) {
            auto id = sg.oldnames_to_ids[n];
            std::set<nodeVisitor> results;
            results.emplace(id, 0, 0);
            results.emplace(-id, 0, 0);
            do {
                if (sglib::OutputLogLevel >= sglib::DEBUG) {
                    sglib::OutputLog(sglib::DEBUG) << "To visit nodes ";
                    for (const auto &n:results) {
                        sglib::OutputLog(sglib::DEBUG, false) << n.node << ":" << n.path_length << ", ";
                    }
                    sglib::OutputLog(sglib::DEBUG, false) << std::endl;
                }
                auto toVisit = *results.begin();
                results.erase(toVisit);
                // Explore node forward
                sglib::OutputLog(sglib::DEBUG) << "Visiting " << toVisit.node << ":" << toVisit.path_length << std::endl;
                auto visitedNodes = sg.depth_first_search(toVisit, size_limit, edge_limit, resultNodes);

                if (sglib::OutputLogLevel >= sglib::DEBUG) {
                    sglib::OutputLog(sglib::DEBUG) << "Visited nodes ";
                    for (const auto &n:visitedNodes) {
                        sglib::OutputLog(sglib::DEBUG, false) << n.node << ":" << n.path_length << ", ";
                    }
                    sglib::OutputLog(sglib::DEBUG, false) << std::endl;
                }

                // For each visited node
                for (const auto &node:visitedNodes) {
                    // If is the same as the node I am visiting, do nothing
                    if (node == toVisit) continue;

                    auto visitedNode(resultNodes.find(node));

                    // If was in previous results
                    if (visitedNode != resultNodes.end()) {
                        // If now has a shorter path or closer distance
                        if (node.path_length > visitedNode->path_length or node.dist > visitedNode->dist) {
                            resultNodes.erase(visitedNode);
                        } else {
                            continue;
                        }
                    }
                    sglib::OutputLog(sglib::DEBUG) << "Inserting node to visit: " << node << std::endl;
                    results.insert(node);
                    nodeVisitor rev = nodeVisitor(-1*node.node, node.dist, node.path_length);
                    sglib::OutputLog(sglib::DEBUG) << "Inserting node to visit: " << rev << std::endl;
                    results.insert(rev);
                    resultNodes.emplace(node);
                    if (sglib::OutputLogLevel >= sglib::DEBUG) {
                        sglib::OutputLog(sglib::DEBUG) << "Result nodes ";
                        for (const auto &n:resultNodes) {
                            sglib::OutputLog(sglib::DEBUG, false) << n.node << ":" << n.path_length << ", ";
                        }
                        sglib::OutputLog(sglib::DEBUG, false) << std::endl;
                    }

                }
                // While exploration results > 0 explore resulting nodes
            } while (!results.empty());

        }

        sglib::OutputLog() << "Nodes in solution" << std::endl;
        for (const auto &n: resultNodes) {
            sglib::OutputLog(sglib::INFO, false) << sg.oldnames[std::abs(n.node)] << " ";
        }
        sglib::OutputLog(sglib::INFO,false) << std::endl;
        sglib::OutputLog() << resultNodes.size() << " nodes in solution\n";
        std::vector<sgNodeID_t > subnodes;
        for (const auto &n:resultNodes) {
            subnodes.emplace_back(n.node);
        }
        SequenceSubGraph ssg(sg, subnodes);
        ssg.write_to_gfa(output_prefix+"subgraph.gfa");
    }

    if (!query_file.empty()) {

    }

    sglib::OutputLog() << "Done" << std::endl;
}