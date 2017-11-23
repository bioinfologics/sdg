//
// Created by Katie Barr (EI) on 13/11/2017.
//

#include <iostream>
#include <fstream>
#include "cxxopts.hpp"
#include "sglib/SequenceGraph.hpp"
#include "sglib/HaplotypeScorer.hpp"
#include "sglib/PhaseScaffolder.h"


int main(int argc, char * argv[]) {
    std::string gfa_filename,bubble_contigs_filename,output_prefix, reads1,reads2;
    std::vector<std::string>  dump_mapped, load_mapped;

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
                ("1,read1", "input reads, left", cxxopts::value<std::string>(reads1))
                ("2,read2", "input reads, right", cxxopts::value<std::string>(reads2))
                ("d,dump_to", "dump mapped reads to file", cxxopts::value<std::vector<std::string>>(dump_mapped))
                ("l,load_from", "load mapped reads from file", cxxopts::value<std::vector<std::string>>(load_mapped));


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


    std::cout<< "Welcome to phaser"<<std::endl<<std::endl;

    PhaseScaffolder ps = PhaseScaffolder(gfa_filename);
    ps.load_mappings(reads1, reads2);
    ps.phase_components();
    SequenceGraph sg;
    sg.load_from_gfa(gfa_filename);
    std::cout << sg.oldnames_to_ids.size() << " sg oodes size: " << sg.nodes.size() << std::endl;
    std::cout << "Edge 0: " << sg.oldnames_to_ids["edge0"] << " " << sg.oldnames_to_ids["edge0+"] << std::endl;



    /*(for (auto node: sg.nodes){
        std::cout << node.sequence << " " << std::endl;
    }*/
    // currently takes gfa, file of possible phasings, and reads
    // loads graph, loads phasings, maps reads and computes support:
    /**
     * \todo Load entire GFA rather than subcomponent
     * \todo map reads to entire GFA
    * \todo phase components with > 1 bubbles in turn

     * \todo from bubble contigs, compute possible haplotypes - if too many, split so exponential growth doesn't kill us
     * \todo decide heurisitcs for picking barcode winner, overall winner, when to
     * with all phased, intersect barocdes supporting phasings on different contigs to build up phase blocks

 */


    return 0;
}
