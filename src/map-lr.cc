//
// Created by Luis Yanes (EI) on 26/01/2018.
//

#include <iostream>
#include <fstream>
#include <sglib/factories/ContigBlockFactory.h>
#include <sglib/mappers/LongReadMapper.h>
#include "cxxopts.hpp"
#include "sglib/SequenceGraph.h"
#include "sglib/Scaffolder.hpp"
#include "sglib/GraphPartitioner.hpp"


int main(int argc, char * argv[]) {
    std::string gfa_filename, bubble_contigs_filename, output_prefix, long_reads;
    std::string dump_mapped, load_mapped;
    unsigned int log_level(4);
    uint64_t max_mem_gb(4);
    bool stats_only(false);
    try {
//@formatter:off
        cxxopts::Options options("lr_repeat_resolver", "Long read repeat resolution");
        options.add_options()
                ("help", "Print help", cxxopts::value<std::string>(),"")
                ("g,gfa", "input gfa file", cxxopts::value<std::string>(gfa_filename), "filepath")
                ("o,output", "output file prefix", cxxopts::value<std::string>(output_prefix), "path")
                ("log_level", "output log level", cxxopts::value<unsigned int>(log_level), "uint");
        options.add_options("Long Read Options")
                ("r,long_reads", "input long reads", cxxopts::value<std::string>(long_reads), "filepath")
                ("d,dump_to","dump mapped reads to file",cxxopts::value<std::string>(dump_mapped), "filepath")
                ("l,load_from", "load mapped reads from file", cxxopts::value<std::string>(load_mapped), "filepath")
                ("max_mem", "maximum_memory when mapping (GB, default: 4)", cxxopts::value<uint64_t>(max_mem_gb)->default_value("4"), "GB");
//@formatter:on


        auto result = options.parse(argc, argv);

        if (result.count("help")) {
            std::cout << options.help({""}) << std::endl;
            exit(0);
        }

        if (result.count("g") != 1 or result.count("o") != 1 ) {
            throw cxxopts::OptionException(" please specify input files and output prefix");
        }
        if (long_reads.empty()) {
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

    sglib::OutputLog() << "Welcome to map-lr" << std::endl << std::endl;
    if (gfa_filename.size() <= 4 or gfa_filename.substr(gfa_filename.size() - 4, 4) != ".gfa") {

        throw std::invalid_argument("filename of the gfa input does not end in gfa, it ends in '" +
                                    gfa_filename.substr(gfa_filename.size() - 4, 4) + "'");
    }
    sglib::OutputLogLevel = static_cast<sglib::LogLevels>(log_level);
    max_mem_gb *= GB;
    SequenceGraph sg;
    sg.load_from_gfa(gfa_filename);
    std::vector<PairedReadMapper> mappers;
    mappers.emplace_back(sg);
    mappers.emplace_back(sg);

    unsigned int K = 15;
    uint16_t min_matches = 2;

    LongReadMapper rm(K, sg);
    rm.map_reads(min_matches, long_reads);

    auto repeatyNodes (sg.find_canonical_repeats());
    std::vector<SequenceGraphPath> solvable_paths;
    // For each repeaty node

    unsigned int loops = 0;
    for (const auto &central_repeat_node:repeatyNodes) {
        sglib::OutputLog() << "* Evaluating central node " << central_repeat_node << ": " << std::endl;
        // Find AA,AB,BA,BB nodes
        auto forward_links = sg.get_fw_links(central_repeat_node);
        auto backward_links = sg.get_bw_links(central_repeat_node);
        auto fwdA = forward_links[0].dest;
        auto fwdB = forward_links[1].dest;
        auto bwdA = backward_links[0].dest;
        auto bwdB = backward_links[1].dest;

        // Check it's not a loop
        std::array<sgNodeID_t, 4> nodes = {fwdA, fwdB, bwdA, bwdB};
        if (sg.is_loop(nodes)) {
            loops++;
            continue;
        }
        // Find the reads that cross those nodes
        std::set<uint32_t> readIDs;

        // Define the AA,BB or AB,BA sets and validate they don't cross over
        sglib::OutputLog() <<"Repeat: [ "<< -bwdA <<" | "<< -bwdB <<" ] <-> "<< central_repeat_node <<" <-> [ "<< fwdA <<" | " << fwdB <<" ]"<<std::endl;

    }
    sglib::OutputLog() << repeatyNodes.size() << " repeat nodes generated " << solvable_paths.size() << " solvable paths" << std::endl;
    sglib::OutputLog() << loops << " loops found in repeat nodes" << std::endl;
    sg.write_to_gfa("salida.gfa");
}
