//
// Created by Luis Yanes (EI) on 26/01/2018.
//

#include <iostream>
#include <fstream>
#include <sglib/factories/ContigBlockFactory.h>
#include <sglib/mappers/LongReadMapper.h>
#include "cxxopts.hpp"


int main(int argc, char * argv[]) {
    std::string gfa_filename, bubble_contigs_filename, output_prefix, long_reads;
    std::string dump_mapped, load_mapped;
    unsigned int log_level(4);
    uint64_t max_mem_gb(4);
    bool stats_only(false);
    try {
//@formatter:off
        cxxopts::Options options("map-lr", "LongRead Mapper");
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

    unsigned int K = 15;
    uint16_t min_matches = 2;

    LongReadMapper rm(K, 5, sg);
    rm.map_reads(long_reads, 500);
    /*
     * Print the length of the node and the number of mapped reads
     */

    auto repeatyNodes (sg.find_canonical_repeats());
    // For each repeaty node
    unsigned int loops = 0;
    std::vector<SequenceGraphPath> paths_solved;
    for (const auto &central_repeat_node:repeatyNodes) {
        sglib::OutputLog() << "* Evaluating central node " << central_repeat_node << ": " << std::endl;
        // Find AA,AB,BA,BB nodes
        auto forward_links = sg.get_fw_links(central_repeat_node);
        auto backward_links = sg.get_bw_links(central_repeat_node);
        auto nodeIDs = std::array<sgNodeID_t , 4>
                {forward_links[0].dest, forward_links[1].dest, // fwdA,fwdB
                 backward_links[0].dest,backward_links[1].dest}; // bwdA,bwdB
        if (sg.is_loop(nodeIDs)) {
            loops++;
            sglib::OutputLog() << "This repeat is loopy" << std::endl;
            continue;
        }
        // Find the reads that cross those nodes
        std::array<std::set<uint32_t>,4> readIDs(rm.getReadSets(nodeIDs));
        int min_coverage(3);
        int set_distance_ratio(10);
        // Check set AA,BB
        std::set<uint32_t > AA,AB,BA,BB;
        std::set_intersection(readIDs[0].cbegin(),readIDs[0].cend(),readIDs[2].cbegin(),readIDs[2].cend(), std::inserter(AA,AA.end()));
        std::set_intersection(readIDs[0].cbegin(),readIDs[0].cend(),readIDs[3].cbegin(),readIDs[3].cend(), std::inserter(AB,AB.end()));
        std::set_intersection(readIDs[1].cbegin(),readIDs[1].cend(),readIDs[2].cbegin(),readIDs[2].cend(), std::inserter(BA,BA.end()));
        std::set_intersection(readIDs[1].cbegin(),readIDs[1].cend(),readIDs[3].cbegin(),readIDs[3].cend(), std::inserter(BB,BB.end()));

        sglib::OutputLog() <<"Long reads in "
                << "fA(" << sg.nodes[std::abs(nodeIDs[0])].sequence.length() << "): " << ""<<readIDs[0].size() << " "
                << "fB(" << sg.nodes[std::abs(nodeIDs[1])].sequence.length() << "): " << ""<<readIDs[1].size() << " "
                << "bA(" << sg.nodes[std::abs(nodeIDs[2])].sequence.length() << "): " << ""<<readIDs[2].size() << " "
                << "bB(" << sg.nodes[std::abs(nodeIDs[3])].sequence.length() << "): " << ""<<readIDs[3].size() << "\n";
        sglib::OutputLog() <<"Paths support AA: "<<AA.size()<<"  BB: "<<BB.size()<<"  AB: "<<AB.size()<<"  BA: "<<BA.size() << std::endl;

        if (AA.size() > min_coverage and BB.size() > min_coverage and
            std::min(AA.size(), BB.size()) > set_distance_ratio * std::max(AB.size(), BA.size())) {
            sglib::OutputLog() << " Solved as AA BB !!!" << std::endl;
            paths_solved.push_back(SequenceGraphPath(sg, {-nodeIDs[2], central_repeat_node, nodeIDs[0]}));  // AA
            paths_solved.push_back(SequenceGraphPath(sg, {-nodeIDs[3], central_repeat_node, nodeIDs[1]}));  // BB
        } else if (BA.size() > min_coverage and AB.size() > min_coverage and
                   std::min(BA.size(), AB.size()) > set_distance_ratio * std::max(AA.size(), BB.size())) {
            sglib::OutputLog() << " Solved as AB BA !!!" << std::endl;
            paths_solved.push_back(SequenceGraphPath(sg, {-nodeIDs[3], central_repeat_node, nodeIDs[0]}));  // AB
            paths_solved.push_back(SequenceGraphPath(sg, {-nodeIDs[2], central_repeat_node, nodeIDs[1]}));  // BA
        }
        // Define the AA,BB or AB,BA sets and validate they don't cross over
        sglib::OutputLog() <<"Repeat: [ "<< -nodeIDs[2] <<" | "<< -nodeIDs[3] <<" ] <-> "<< central_repeat_node <<" <-> [ "<< nodeIDs[0] <<" | " << nodeIDs[1] <<" ]"<<std::endl;

    }
    sglib::OutputLog() << repeatyNodes.size() << " repeat nodes generated " << paths_solved.size() << " solvable paths" << std::endl;
    sglib::OutputLog() << loops << " loops found in repeat nodes" << std::endl;
    sg.write_to_gfa("salida.gfa");
}
