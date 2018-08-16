#include <iostream>
#include <fstream>
#include <sglib/datastores/LongReadsDatastore.hpp>
#include "cxxopts.hpp"


int main(int argc, char * argv[]) {
    std::cout << "lr-datastore"<<std::endl<<std::endl;
    std::cout << "Git origin: " << GIT_ORIGIN_URL << " -> "  << GIT_BRANCH << std::endl;
    std::cout << "Git commit: " << GIT_COMMIT_HASH << std::endl<<std::endl;
    std::cout << "Executed command:"<<std::endl;
    for (auto i=0;i<argc;i++) std::cout<<argv[i]<<" ";
    std::cout<<std::endl<<std::endl;


    std::string reads_filename;
    std::string read_type;
    std::string output;

    cxxopts::Options options("lr-datastore", "make lr-datastore");
    options.add_options()
            ("help", "Print help", cxxopts::value<std::string>()->implicit_value(""),"uint")
            ("r,read", "input reads", cxxopts::value<std::string>(reads_filename), "filename")
            ("t,read_type", "One of: pacbio,nanopore", cxxopts::value<std::string>(read_type)->default_value("nanopore"), "pacbio, nanopore")
            ("o,output", "output file", cxxopts::value<std::string>(output), "filename");
    try {
        auto result=options.parse(argc,argv);
        if (result.count("help")) {
            std::cout << options.help({""}) << std::endl;
            exit(0);
        }

        if (reads_filename== "" or read_type == "" or output == "") {
            throw cxxopts::OptionException(" please specify input files, type and output prefix");
        }


    } catch (const cxxopts::OptionException &e) {
        std::cout << "Error parsing options: " << e.what() << std::endl << std::endl
                  << "Use option --help to check command line arguments." << std::endl;
        exit(1);
    }

    LongReadsDatastore lr_ds(reads_filename, output);

    exit(0);
}

