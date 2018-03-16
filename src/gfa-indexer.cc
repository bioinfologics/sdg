#include <iostream>
#include <fstream>
#include <sglib/SequenceGraph.h>
#include <cxxopts.hpp>
#include <sglib/logger/OutputLog.h>
#include <sglib/indexers/skipMerIndex.hpp>


int main(int argc, char * argv[]) {
    std::cout << "Welcome to gfa-indexer"<<std::endl<<std::endl;
    std::cout << "Git origin: " << GIT_ORIGIN_URL << " -> "  << GIT_BRANCH << std::endl;
    std::cout << "Git commit: " << GIT_COMMIT_HASH << std::endl<<std::endl;
    std::cout << "Executed command:"<<std::endl;
    for (auto i=0;i<argc;i++) std::cout<<argv[i]<<" ";
    std::cout<<std::endl<<std::endl;

    std::string gfa_filename,output_prefix;
    uint64_t max_mem_gb=4;
    sglib::OutputLogLevel=sglib::LogLevels::DEBUG;
    cxxopts::Options options("gfa-indexer", "GFA Indexing tool");

    options.add_options()
            ("help", "Print help")
            ("g,gfa", "input gfa file", cxxopts::value<std::string>(gfa_filename))
            ("o,output", "output file prefix", cxxopts::value<std::string>(output_prefix));
    try
    {
        auto result(options.parse(argc, argv));

        if (result.count("help"))
        {
            std::cout << options.help({""}) << std::endl;
            exit(0);
        }

        if (result.count("g")!=1 or result.count("o")!=1) {
            throw cxxopts::OptionException(" please specify input files and output prefix");
        }
    } catch (const cxxopts::OptionException& e)
    {
        std::cerr << "Error parsing options: " << e.what() << std::endl << std::endl;
        std::cout << options.help({""}) << std::endl;
        exit(1);
    }

    std::cout<<std::endl<<"=== Loading GFA ==="<<std::endl;
    SequenceGraph sg;
    sg.load_from_gfa(gfa_filename);

    SkipMerIndex skmNDX(sg, 15, 1, 3, 100);

    return 0;
}

