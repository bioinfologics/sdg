#include <iostream>
#include <fstream>
#include "cxxopts.hpp"
#include "sglib/SequenceGraph.h"


int main(int argc, char * argv[]) {
    std::string gfa_filename,output_prefix,help;
    std::vector<std::string> paths;
    bool stats_only=0;

    try
    {
        cxxopts::Options options("gfa-dumppath", "GFA dump path tool");

        options.add_options()
                ("help", "Print help", cxxopts::value<std::string>(help), "")
                ("g,gfa", "input gfa file", cxxopts::value<std::string>(gfa_filename), "GFA")
                ("p,path", "coma separated lists of nodes", cxxopts::value<std::vector<std::string>>(paths), "ex: 123-,34+,16+")
                ("o,output", "output file prefix", cxxopts::value<std::string>(output_prefix), "path");



        auto result = options.parse(argc, argv);

        if (result.count("help"))
        {
            std::cout << options.help({""}) << std::endl;
            exit(0);
        }

        if (result.count("g")!=1 or result.count("p")<1 or result.count("o")!=1) {
            throw cxxopts::OptionException(" please specify input files and output prefix");
        }



    } catch (const cxxopts::OptionException& e)
    {
        std::cout << "Error parsing options: " << e.what() << std::endl << std::endl
                <<"Use option --help to check command line arguments." << std::endl;
        exit(1);
    }


    std::cout<< "Welcome to gfaqc"<<std::endl<<std::endl;

    SequenceGraph sg;
    sg.load_from_gfa(gfa_filename);
    std::ofstream ofasta(output_prefix+".fa");
    for (auto &p:paths){
        std::cout<<"Extracting path: "<<p<<std::endl;
        try {
            SequenceGraphPath path(sg,sg.oldnames_to_nodes(p));
            std::cout<<"Path in current Ids: ";
            for (auto n:path.nodes) {std::cout<<" "<<n;}
            std::cout<<std::endl;
            ofasta<<path.get_fasta_header()<<std::endl;
            ofasta<<path.get_sequence()<<std::endl;

        } catch (const std::exception& e) {
            std::cout << "Can't dump path "<<p<<", check path is valid and files are accessible"<<std::endl;
        }
    }
    return 0;
}

