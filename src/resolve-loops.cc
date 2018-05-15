//
// Created by Luis Yanes (EI) on 26/01/2018.
//

#include <iostream>
#include <fstream>
#include <vector>
#include <cxxopts.hpp>
#include <sglib/logger/OutputLog.h>
#include <sglib/filesystem/helpers.h>
#include <sglib/graph/SequenceGraph.hpp>


int main(int argc, char * argv[]) {
    sglib::OutputLog(false) << "Welcome to Loop resolution"<<std::endl<<std::endl;
    sglib::OutputLog(false) << "Git origin: " << GIT_ORIGIN_URL << " -> "  << GIT_BRANCH << std::endl;
    sglib::OutputLog(false) << "Git commit: " << GIT_COMMIT_HASH << std::endl<<std::endl;
    sglib::OutputLog() << "Executed command:"<<std::endl;
    for (auto i=0;i<argc;i++) sglib::OutputLog(false) <<argv[i]<<" ";
    sglib::OutputLog(false) <<std::endl<<std::endl;

    std::string gfa_filename, bubble_contigs_filename, output_prefix, lr_datastore;
    std::string dump_mapped, load_mapped;
    unsigned int log_level(0);
    uint64_t max_mem_gb(4);
    bool stats_only(false);
    uint8_t K(15), W(0);
//@formatter:off
    cxxopts::Options options("resolve-loops", "Loop resolution");
    options.add_options()
            ("help", "Print help", cxxopts::value<std::string>(),"")
            ("g,gfa", "input gfa file", cxxopts::value<std::string>(gfa_filename), "filepath")
            ("o,output", "output file prefix", cxxopts::value<std::string>(output_prefix), "path")
            ("l,log_level", "output log level", cxxopts::value<unsigned int>(log_level)->default_value("0"), "uint")
            ("d,datastore", "long read datastore", cxxopts::value<std::string>(lr_datastore), "filepath");
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
        if (gfa_filename.size() <= 4 or gfa_filename.substr(gfa_filename.size() - 4, 4) != ".gfa") {

            throw cxxopts::OptionException("filename of the gfa input does not end in gfa, it ends in '" +
                                        gfa_filename.substr(gfa_filename.size() - 4, 4) + "'");
        }
    } catch (const cxxopts::OptionException &e) {
        std::cout << "Error parsing options: " << e.what() << std::endl;
        std::cout << options.help({""}) << std::endl;
        exit(1);
    }
    if (!sglib::check_or_create_directory(output_prefix)) {
        exit(1);
    }

    sglib::OutputLogLevel = static_cast<sglib::LogLevels>(log_level);
    SequenceGraph sg;
    sg.load_from_gfa(gfa_filename);

    unsigned int complexity=2;
    // If exploring complexity nodes I return to the same node, then it is a loop.
    const auto loopy_nodes(sg.get_loopy_nodes(complexity));
    unsigned int loop(0);
    unsigned int validLoops(0);
    unsigned int smallFlanks(0);
    unsigned int missingFlank(0);
    std::vector<sgNodeID_t > validLoopNodes;
    validLoopNodes.reserve(1000);
    for (const auto &loopy_node:loopy_nodes) {
        bool isValidLoop(true);
        loop++;
        const auto flank(sg.get_flanking_nodes(loopy_node));
        if (flank.size() <= 1) {isValidLoop=false; missingFlank++;}
        if (!isValidLoop) continue;
        for (const auto &fn: flank ) {
            if (sg.nodes[(fn>0?fn:-fn)].sequence.size() < 1000) {
                isValidLoop = false;
                smallFlanks++;
                break;
            }
        }
        if (!isValidLoop) continue;
        validLoopNodes.emplace_back(loopy_node);
        std::cout << "Resolving loop " << loop << " / " << loopy_nodes.size() << "\n";
        std::cout << loopy_node << "\t" << sg.oldnames[loopy_node] << "\n";
        std::cout << "Flanks = ";
        std::copy(flank.begin(),flank.end(),std::ostream_iterator<sgNodeID_t>(std::cout, ", "));
        std::cout << std::endl;
    }
    sglib::OutputLog() << "Loops missing flanks = " << missingFlank << std::endl;
    sglib::OutputLog() << "Flanks smaller than 1000bp = " << smallFlanks << std::endl;

    return 0;
}
