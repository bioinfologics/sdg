//
// Created by Luis Yanes (EI) on 25/01/2018.
//

#include <iostream>
#include <fstream>
#include "cxxopts.hpp"
#include "sglib/SequenceGraph.h"
#include "sglib/HaplotypeScorer.hpp"
#include "sglib/PhaseScaffolder.h"


int main(int argc, char * argv[]) {
    std::string gfa_filename,output_prefix,reads1,reads2;
    std::vector<std::string>  dump_mapped, load_mapped;
    uint64_t max_mem_gb=4;
    bool stats_only=0;

    try
    {
        cxxopts::Options options("repeat_resolver", "Tool to expand canonical repeats");

        options.add_options()
                ("help", "Print help")
                ("g,gfa", "input gfa file", cxxopts::value<std::string>(gfa_filename))
                ("o,output", "output file prefix", cxxopts::value<std::string>(output_prefix));
        options.add_options("10x reads options")
                ("1,read1", "input reads, left", cxxopts::value<std::string>(reads1))
                ("2,read2", "input reads, right", cxxopts::value<std::string>(reads2))
                ("d,dump_to", "dump mapped reads to file", cxxopts::value<std::vector<std::string>>(dump_mapped))
                ("l,load_from", "load mapped reads from file", cxxopts::value<std::vector<std::string>>(load_mapped))
                ("max_mem", "maximum_memory when mapping (GB, default: 4)", cxxopts::value<uint64_t>(max_mem_gb));

        auto result = options.parse(argc, argv);

        if (result.count("help"))
        {
            std::cout << options.help({""}) << std::endl;
            exit(0);
        }

        if (result.count("g")!=1 or result.count("o")!=1 or result.count("c") != 1) {
            throw cxxopts::OptionException(" please specify input files and output prefix");
        }
        if (result.count("1")!=1 or result.count("2")!=1) {
            throw cxxopts::OptionException(" please specify  R1 and R2 of 10x readfiles ");

        }

    } catch (const cxxopts::OptionException& e)
    {
        std::cout << "Error parsing options: " << e.what() << std::endl << std::endl
                  <<"Use option --help to check command line arguments." << std::endl;
        exit(1);
    }
    if (!sglib::check_or_create_directory(output_prefix)) {
        exit(1);
    }

    std::cout<< "Welcome to repeat_resolver"<<std::endl<<std::endl;
    if (gfa_filename.size()<=4 or gfa_filename.substr(gfa_filename.size()-4,4)!=".gfa") {

        throw std::invalid_argument("filename of the gfa input does not end in gfa, it ends in '" +
                                    gfa_filename.substr(gfa_filename.size() - 4, 4) + "'");
    }
    max_mem_gb *= GB;
    SequenceGraph sg;

    sg.load_from_gfa(gfa_filename);
    PairedReadMapper mapper(sg);

    std::ifstream mappings(output_prefix+"mappings.map");
    if (mappings.good()) {
        mapper.load_from_disk(output_prefix+"mappings.map");
    } else {
        mapper.map_reads(reads1, reads2, PairedReadMapper::prm10x, max_mem_gb);
        mapper.save_to_disk(output_prefix + "mappings.map");
    }

    return 0;
}
