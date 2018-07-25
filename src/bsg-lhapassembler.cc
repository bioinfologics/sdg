#include <iostream>
#include <fstream>
#include <sglib/WorkSpace.hpp>
#include <sglib/processors/LocalHaplotypeAssembler.hpp>
#include "sglib/logger/OutputLog.h"
#include "cxxopts.hpp"


int main(int argc, char * argv[]) {
    std::cout << "Welcome to bsg-lhapassembler"<<std::endl<<std::endl;
    std::cout << "Git origin: " << GIT_ORIGIN_URL << " -> "  << GIT_BRANCH << std::endl;
    std::cout << "Git commit: " << GIT_COMMIT_HASH << std::endl<<std::endl;
    std::cout << "Executed command:"<<std::endl;
    for (auto i=0;i<argc;i++) std::cout<<argv[i]<<" ";
    std::cout<<std::endl<<std::endl;

    sglib::OutputLogLevel=sglib::LogLevels::DEBUG;

    std::string workspace_file,problem_file,output_prefix;

    uint8_t k=63;
    int min_cvg=5;
    try
    {
        cxxopts::Options options("bsg-lhapassembler", "Local Haplotype(specific) Assembler");

        options.add_options()
                ("help", "Print help")
                ("w,workspace", "input workspace", cxxopts::value<std::string>(workspace_file))
                ("p,problem", "problem file", cxxopts::value<std::string>(problem_file))
                ("o,output", "output file prefix", cxxopts::value<std::string>(output_prefix));

        options.add_options("Assembly options")
                ("k,kmer_size","k for local assembly",cxxopts::value<uint8_t>(k))
                ("c,min_cvg","minimum coverga for local assemblies",cxxopts::value<int>(min_cvg));


        auto result(options.parse(argc, argv));

        if (result.count("help"))
        {
            std::cout << options.help({"","Backbone (unique anchors linkage)","Untangling functions","Development"}) << std::endl;
            exit(0);
        }

        if (result.count("p")!=1 or result.count("o")!=1) {
            throw cxxopts::OptionException(" please specify input problem and output prefix");
        }



    } catch (const cxxopts::OptionException& e)
    {
        std::cout << "Error parsing options: " << e.what() << std::endl << std::endl
                <<"Use option --help to check command line arguments." << std::endl;
        exit(1);
    }

    std::cout<<std::endl;


    //======= WORKSPACE LOAD AND CHECKS ======
    WorkSpace ws;
    LocalHaplotypeAssembler lha(ws);

    if (!workspace_file.empty()) {
        sglib::OutputLog() << "Loading Workspace..." << std::endl;
        ws.load_from_disk(workspace_file);
        if (!ws.sg.is_sane()) {
            sglib::OutputLog() << "ERROR: sg.is_sane() = false" << std::endl;
            return 1;
        }

        //TODO: other checks? reads mapped to valid nodes and such?

        sglib::OutputLog() << "Loading Workspace DONE" << std::endl;
        lha.init_from_file(problem_file);
    }
    else {
        lha.init_from_full_file(problem_file);
    }


    lha.assemble(k,min_cvg,false,true,output_prefix);

    return 0;
}

