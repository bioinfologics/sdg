#include <iostream>
#include <fstream>
#include "cxxopts.hpp"
#include "sglib/SequenceGraph.h"


int main(int argc, char * argv[]) {
    std::string gfa_filename,ref_filename,output_prefix;
    bool stats_only=0;

    try
    {
        cxxopts::Options options("gfaqc", "GFA QC tool");

        options.add_options()
                ("help", "Print help")
                ("g,gfa", "input gfa file", cxxopts::value<std::string>(gfa_filename))
                ("R,reference", "reference file (FASTA)", cxxopts::value<std::string>(ref_filename))
                ("o,output", "output file prefix", cxxopts::value<std::string>(output_prefix))
                ("s,stats_only", "do not dump detailed information (default=0)", cxxopts::value<bool>(stats_only));



        auto result(options.parse(argc, argv));

        if (result.count("help"))
        {
            std::cout << options.help({""}) << std::endl;
            exit(0);
        }

        if (result.count("g")!=1 or result.count("R")!=1 or result.count("o")!=1) {
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
    sg.write_to_gfa(output_prefix+".gfa");
    return 0;
}

