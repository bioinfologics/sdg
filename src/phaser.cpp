//
// Created by Katie Barr (EI) on 13/11/2017.
//

#include <iostream>
#include <fstream>
#include "deps/cxxopts.hpp"
#include "sglib/SequenceGraph.hpp"
#include "haplotype_scorer.h"


int main(int argc, char * argv[]) {
    std::string gfa_filename,bubble_contigs_filename,output_prefix;
    std::vector<std::string> reads1,reads2, dump_mapped, load_mapped;

    bool stats_only=0;

    try
    {
        cxxopts::Options options("gfaqc", "GFA QC tool");

        options.add_options()
                ("help", "Print help")
                ("g,gfa", "input gfa file", cxxopts::value<std::string>(gfa_filename))
                ("c,contigs", "Bubble contigs to phase", cxxopts::value<std::string>(bubble_contigs_filename))
                ("o,output", "output file prefix", cxxopts::value<std::string>(output_prefix));
        options.add_options("Paired reads options")
                ("1,read1", "input reads, left", cxxopts::value<std::vector<std::string>>(reads1))
                ("2,read2", "input reads, right", cxxopts::value<std::vector<std::string>>(reads2))
                ("d,dump_to", "dump mapped reads to file", cxxopts::value<std::vector<std::string>>(dump_mapped))
                ("l,load_from", "load mapped reads from file", cxxopts::value<std::vector<std::string>>(load_mapped));


        options.parse(argc, argv);

        if (options.count("help"))
        {
            std::cout << options.help({""}) << std::endl;
            exit(0);
        }

        if (options.count("g")!=1 or options.count("o")!=1 or options.count("c") != 1) {
            throw cxxopts::OptionException(" please specify input files and output prefix");
        }
        if (options.count("1")!=1 or options.count("2")!=1) {
            throw cxxopts::OptionException(" please specify  R1 and R2 of 10x readfiles ");

        }




    } catch (const cxxopts::OptionException& e)
    {
        std::cout << "Error parsing options: " << e.what() << std::endl << std::endl
                  <<"Use option --help to check command line arguments." << std::endl;
        exit(1);
    }


    std::cout<< "Welcome to phaser"<<std::endl<<std::endl;

    SequenceGraph sg;
    sg.load_from_gfa(gfa_filename);
    std::cout << sg.oldnames_to_ids.size() << std::endl;
    std::cout << "Edge 0: " << sg.oldnames_to_ids["edge0"] << " " << sg.oldnames_to_ids["edge0+"] << std::endl;
    // take list of contigs to phase from input
    // generate possible haplotypes
    // load mappings from input or disk
    // score haplotypes
    // -----------------------------
    // instead of loading gfa and contig list separately, start by using gfa of subcomponent as input
    // be given names of bubble contigs to avoid having to code bubble finding now
    HaplotypeScorer hs(sg, "blah.txt");
    // for now supply next file with list of possible phasings
    hs.load_haplotypes(bubble_contigs_filename, 2);
    sg.write_to_gfa(output_prefix+".gfa");
    return 0;
}
