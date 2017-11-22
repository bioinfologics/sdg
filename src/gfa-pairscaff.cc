#include <iostream>
#include <fstream>
#include <sglib/PairedReadMapper.hpp>
#include <sglib/Scaffolder.hpp>
#include "sglib/SequenceGraph.hpp"
#include "cxxopts.hpp"


int main(int argc, char * argv[]) {
    std::string gfa_filename,output_prefix;
    std::vector<std::string> reads1,reads2, dump_mapped, load_mapped;
    bool stats_only=0;

    try
    {
        cxxopts::Options options("gfa-pairscaff", "GFA paired reads scaffolder");

        options.add_options()
                ("help", "Print help")
                ("g,gfa", "input gfa file", cxxopts::value<std::string>(gfa_filename))
                ("o,output", "output file prefix", cxxopts::value<std::string>(output_prefix));
        options.add_options("Paired reads options")
                ("1,read1", "input reads, left", cxxopts::value<std::vector<std::string>>(reads1))
                ("2,read2", "input reads, right", cxxopts::value<std::vector<std::string>>(reads2))
                ("d,dump_to", "dump mapped reads to file", cxxopts::value<std::vector<std::string>>(dump_mapped))
                ("l,load_from", "load mapped reads from file", cxxopts::value<std::vector<std::string>>(load_mapped));


        auto result(options.parse(argc, argv));

        if (result.count("help"))
        {
            std::cout << options.help({"","Paired reads options"}) << std::endl;
            exit(0);
        }

        if (result.count("g")!=1 or (result.count("1")<1 and result.count("l")<1) or result.count("o")!=1) {
            throw cxxopts::OptionException(" please specify input files and output prefix");
        }

        if ( result.count("1")!=result.count("2")){
            throw cxxopts::OptionException(" please specify read1 and read2 files in pairs");
        }

        if ( result.count("d")>0 and result.count("d")!=result.count("2")){
            throw cxxopts::OptionException(" please specify dump files for all libs, or for none");
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
    {
        std::ofstream ccf(output_prefix + "_initial_components.txt");
        uint16_t i = 1;
        std::map<sgNodeID_t, std::string> nodes_to_oldnames;
        for (auto &nm:sg.oldnames_to_ids) nodes_to_oldnames[(nm.second > 0 ? nm.second : -nm.second)] = nm.first;
        for (auto &c:sg.connected_components()) {
            ccf << "Component #" << i++ << ":";
            for (auto &n:c) ccf << " " << nodes_to_oldnames[n];
            ccf << std::endl;
        }
    }


    //Now try read mapping (as of now, just the first library)

    std::vector<PairedReadMapper> mappers;
    for(auto loadfile:load_mapped){
        mappers.emplace_back(sg);
        mappers.back().load_from_disk(loadfile);
        mappers.back().print_stats();
    }
    for(int lib=0;lib<reads1.size();lib++) {
        mappers.emplace_back(sg);
        mappers.back().map_reads(reads1[lib], reads2[lib], prmPE);
        mappers.back().print_stats();
        if (dump_mapped.size() > 0) {
            std::cout<<"dumping map to "<<dump_mapped[lib]<<std::endl;
            mappers.back().save_to_disk(dump_mapped[lib]);
        }
    }
    //simple scaffolder:

    Scaffolder scaff(sg,mappers);
    scaff.find_canonical_repeats();






    sg.write_to_gfa(output_prefix+"_scaffolded.gfa");
    return 0;
}

