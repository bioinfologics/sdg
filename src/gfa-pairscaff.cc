#include <iostream>
#include <fstream>
#include <sglib/PairedReadMapper.hpp>
#include <sglib/Scaffolder.hpp>
#include <sglib/KmerCompressionIndex.hpp>
#include "sglib/SequenceGraph.hpp"
#include "cxxopts.hpp"


int main(int argc, char * argv[]) {
    std::cout << "Welcome to gfa-pairscaff"<<std::endl<<std::endl;
    std::cout << "Git origin: " << GIT_ORIGIN_URL << " -> "  << GIT_BRANCH << std::endl;
    std::cout << "Git commit: " << GIT_COMMIT_HASH << std::endl<<std::endl;

    std::string gfa_filename,output_prefix, load_cidx, dump_cidx;
    std::vector<std::string> reads1,reads2,reads_type,cidxreads1,cidxreads2, dump_mapped, load_mapped;
    bool stats_only=0;
    uint64_t max_mem_gb=4;

    try
    {
        cxxopts::Options options("gfa-pairscaff", "GFA paired reads scaffolder");

        options.add_options()
                ("help", "Print help")
                ("g,gfa", "input gfa file", cxxopts::value<std::string>(gfa_filename))
                ("o,output", "output file prefix", cxxopts::value<std::string>(output_prefix));
        options.add_options("Compression Index Options")
                ("cidxread1", "compression index input reads, left", cxxopts::value<std::vector<std::string>>(cidxreads1))
                ("cidxread2", "compression index input reads, right", cxxopts::value<std::vector<std::string>>(cidxreads2))
                ("load_cidx", "load compression index filename", cxxopts::value<std::string>(load_cidx))
                ("dump_cidx", "dump compression index filename", cxxopts::value<std::string>(dump_cidx))
                ;
        options.add_options("Paired reads options")
                ("1,read1", "input reads, left", cxxopts::value<std::vector<std::string>>(reads1))
                ("2,read2", "input reads, right", cxxopts::value<std::vector<std::string>>(reads2))
                ("read_type", "One of: pe,10x", cxxopts::value<std::vector<std::string>>(reads_type))
                ("d,dump_to", "dump mapped reads to file", cxxopts::value<std::vector<std::string>>(dump_mapped))
                ("l,load_from", "load mapped reads from file", cxxopts::value<std::vector<std::string>>(load_mapped))
                ("max_mem", "maximum_memory when mapping (GB, default: 4)", cxxopts::value<uint64_t>(max_mem_gb));


        auto result(options.parse(argc, argv));

        if (result.count("help"))
        {
            std::cout << options.help({"","Paired reads options","Compression Index Options"}) << std::endl;
            exit(0);
        }

        if (result.count("g")!=1 or (result.count("1")<1 and result.count("l")<1) or result.count("o")!=1) {
            throw cxxopts::OptionException(" please specify input files and output prefix");
        }

        if ( result.count("1")!=result.count("2") or result.count("1")!=result.count("read_type")){
            throw cxxopts::OptionException(" please specify read1, read2 and read_type files in triplets");
        }

        if ( result.count("cidxread1")!=result.count("cidxread2")){
            throw cxxopts::OptionException(" please specify cidxread1 and cidxread2 files in pairs");
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


    std::cout << "Executed command:"<<std::endl;
    for (auto i=0;i<argc;i++) std::cout<<argv[i]<<" ";
    std::cout<<std::endl<<std::endl;




    std::cout<<std::endl<<"=== Loading GFA ==="<<std::endl;
    if (gfa_filename.size()<=4 or gfa_filename.substr(gfa_filename.size()-4,4)!=".gfa") {

        throw std::invalid_argument("filename of the gfa input does not end in gfa, it ends in '" +
                                    gfa_filename.substr(gfa_filename.size() - 4, 4) + "'");
    }
    auto fasta_filename=gfa_filename.substr(0,gfa_filename.size()-4)+".fasta";
    SequenceGraph sg;
    sg.load_from_gfa(gfa_filename);

    std::cout<<std::endl<<"=== Loading reads compression index ==="<<std::endl;
    //compression index
    KmerCompressionIndex kci(sg,max_mem_gb*1024L*1024L*1024L);
    if (load_cidx!=""){
        kci.load_from_disk(load_cidx);
    } else {
        kci.index_graph();
        for(int lib=0;lib<cidxreads1.size();lib++) {
            kci.start_new_count();
            kci.add_counts_from_file(cidxreads1[lib]);
            kci.add_counts_from_file(cidxreads2[lib]);
        }
    }
    if (dump_cidx!=""){
        kci.save_to_disk(dump_cidx);
    }

    if (kci.read_counts.size()>0) {
        kci.compute_compression_stats();
        kci.dump_histogram("kci_histogram.csv");
    }

    std::cout<<std::endl<<"=== Mapping reads ==="<<std::endl;
    //read mapping/loading
    std::vector<PairedReadMapper> mappers;
    for(auto loadfile:load_mapped){
        mappers.emplace_back(sg);
        mappers.back().load_from_disk(loadfile);
        mappers.back().print_stats();
    }
    for(int lib=0;lib<reads1.size();lib++) {
        mappers.emplace_back(sg);
        if (reads_type[lib]=="10x") {
            mappers.back().map_reads(reads1[lib], reads2[lib], prm10x, max_mem_gb * 1024L * 1024L * 1024L);
        } else {
            mappers.back().map_reads(reads1[lib], reads2[lib], prmPE, max_mem_gb * 1024L * 1024L * 1024L);
        }
        mappers.back().print_stats();
        if (dump_mapped.size() > 0) {
            std::cout<<"dumping map to "<<dump_mapped[lib]<<std::endl;
            mappers.back().save_to_disk(dump_mapped[lib]);
        }
    }

    std::cout<<std::endl<<"=== Scaffolding ==="<<std::endl;

    Scaffolder scaff(sg,mappers,kci);

    scaff.pop_unsupported_shortbubbles();
    sg.join_all_unitigs();
    for (auto &m:mappers) {
        std::cout<<"removing obsolete mappings from "<<m.read1filename<<" and "<<m.read2filename<<std::endl;
        m.remove_obsolete_mappings();
        m.print_stats();
        m.remap_reads();
        m.print_stats();
    }
    //TODO: a lot of repeats are small repeats creating a big "loop", account for those!
    //scaff.find_canonical_repeats();






    sg.write_to_gfa(output_prefix+"_scaffolded.gfa");
    return 0;
}

