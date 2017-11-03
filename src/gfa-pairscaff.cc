#include <iostream>
#include <fstream>
#include "deps/cxxopts.hpp"
#include "sglib/SequenceGraph.hpp"


int main(int argc, char * argv[]) {
    std::string gfa_filename,output_prefix;
    std::vector<std::string> reads1,reads2;
    bool stats_only=0;

    try
    {
        cxxopts::Options options("gfa-pairscaff", "GFA paired reads scaffolder");

        options.add_options()
                ("help", "Print help")
                ("g,gfa", "input gfa file", cxxopts::value<std::string>(gfa_filename))
                ("1,read1", "input reads, left", cxxopts::value<std::vector<std::string>>(reads1))
                ("2,read2", "input reads, right", cxxopts::value<std::vector<std::string>>(reads2))
                ("o,output", "output file prefix", cxxopts::value<std::string>(output_prefix));



        options.parse(argc, argv);

        if (options.count("help"))
        {
            std::cout << options.help({""}) << std::endl;
            exit(0);
        }

        if (options.count("g")!=1 or options.count("1")<1 or options.count("o")!=1) {
            throw cxxopts::OptionException(" please specify input files and output prefix");
        }

        if ( options.count("1")!=options.count("2")){
            throw cxxopts::OptionException(" please specify read1 and read2 files in pairs");
        }



    } catch (const cxxopts::OptionException& e)
    {
        std::cout << "Error parsing options: " << e.what() << std::endl << std::endl
                <<"Use option --help to check command line arguments." << std::endl;
        exit(1);
    }


    std::cout<< "Welcome to gfa-pairscaff"<<std::endl<<std::endl;

    SequenceGraph sg;
    sg.load_from_gfa(gfa_filename);

    //First, report on connected components.



    sg.write_to_gfa(output_prefix+"_scaffolded.gfa");
    return 0;
}

