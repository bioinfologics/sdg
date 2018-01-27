//
// Created by Luis Yanes (EI) on 26/01/2018.
//

#include <iostream>
#include <fstream>
#include "cxxopts.hpp"
#include "sglib/SequenceGraph.hpp"
#include "sglib/HaplotypeScorer.hpp"
#include "sglib/PhaseScaffolder.h"


int main(int argc, char * argv[]) {
    std::string gfa_filename, bubble_contigs_filename, output_prefix, long_reads;
    std::vector<std::string> dump_mapped, load_mapped;
    uint64_t max_mem_gb = 4;
    bool stats_only = 0;

    try {
//@formatter:off
        cxxopts::Options options("map-lr", "Long read mapping tool");
        options.add_options()
                ("help", "Print help", cxxopts::value<std::string>(),"")
                ("g,gfa", "input gfa file", cxxopts::value<std::string>(gfa_filename), "filepath")
                ("o,output", "output file prefix", cxxopts::value<std::string>(output_prefix), "path");
        options.add_options("Read options")
                ("r,long_reads", "input long reads", cxxopts::value<std::string>(long_reads), "filepath")
                ("d,dump_to","dump mapped reads to file",cxxopts::value<std::vector<std::string>>(dump_mapped), "filepath")
                ("l,load_from", "load mapped reads from file", cxxopts::value<std::vector<std::string>>(load_mapped), "filepath")
                ("max_mem", "maximum_memory when mapping (GB, default: 4)", cxxopts::value<uint64_t>(max_mem_gb), "GB");
//@formatter:on


        auto result = options.parse(argc, argv);

        if (result.count("help")) {
            std::cout << options.help({""}) << std::endl;
            exit(0);
        }

        if (result.count("g") != 1 or result.count("o") != 1 or result.count("c") != 1) {
            throw cxxopts::OptionException(" please specify input files and output prefix");
        }
        if (result.count("r") != 1 ) {
            throw cxxopts::OptionException(" please specify a long reads file");

        }


    } catch (const cxxopts::OptionException &e) {
        std::cout << "Error parsing options: " << e.what() << std::endl << std::endl
                  << "Use option --help to check command line arguments." << std::endl;
        exit(1);
    }
    if (!sglib::check_or_create_directory(output_prefix)) {
        exit(1);
    }

    std::cout << "Welcome to map-lr" << std::endl << std::endl;
    if (gfa_filename.size() <= 4 or gfa_filename.substr(gfa_filename.size() - 4, 4) != ".gfa") {

        throw std::invalid_argument("filename of the gfa input does not end in gfa, it ends in '" +
                                    gfa_filename.substr(gfa_filename.size() - 4, 4) + "'");
    }
    max_mem_gb *= GB;
    SequenceGraph sg;
    sg.load_from_gfa(gfa_filename);
    PairedReadMapper mapper(sg);
    mapper.map_reads(long_reads, max_mem_gb);


}