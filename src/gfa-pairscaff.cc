#include <iostream>
#include <fstream>
#include <sglib/PairedReadMapper.h>
#include <sglib/Scaffolder.hpp>
#include <sglib/KmerCompressionIndex.hpp>
#include <sglib/GraphPartitioner.hpp>
#include "sglib/SequenceGraph.h"
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
            kci.add_counts_from_file({cidxreads1[lib]});
            kci.add_counts_from_file({cidxreads2[lib]});
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
            mappers.back().map_reads(reads1[lib], reads2[lib], PairedReadMapper::prm10x, max_mem_gb * 1024L * 1024L * 1024L);
        } else {
            mappers.back().map_reads(reads1[lib], reads2[lib], PairedReadMapper::prmPE, max_mem_gb * 1024L * 1024L * 1024L);
        }
        mappers.back().print_stats();
        if (dump_mapped.size() > 0) {
            std::cout<<"dumping map to "<<dump_mapped[lib]<<std::endl;
            mappers.back().save_to_disk(dump_mapped[lib]);
        }
    }

    std::cout<<std::endl<<"=== Scaffolding ==="<<std::endl;

    std::vector<LinkedReadMapper> linkedReadMapper;
    Scaffolder scaff(sg,mappers, linkedReadMapper, kci);

    std::cout<<std::endl<<"Testing GraphPartitioner"<<std::endl;

    auto bubblies=scaff.get_all_bubbly_subgraphs();
    std::cout<<"Starting with "<<bubblies.size()<<" possible bubbles"<<std::endl;
    uint64_t solved_count=0;
    for (auto bubbly:bubblies) {
        std::cout<<std::endl<<"=== Analysing subgraph with "<<bubbly.nodes.size()<<" nodes"<<std::endl;
        GraphPartitioner partitioner(sg,scaff.rmappers,scaff.kci);
        auto tp=partitioner.tags_patterns(bubbly);
        /*std::cout<<"Tag patterns (from "<<tp.size()<<" tags):"<<std::endl;
        for (auto stp:tp) {
            for (auto p:stp) std::cout<<" "<<(p ? 1:0);
            std::cout<<std::endl;
        }*/
        std::cout<<std::endl;
        auto parts=partitioner.generate_partitions(bubbly,tp);
        std::cout<<"Partitions:"<<std::endl;
        unsigned pnumb=0;
        for (auto &psg:parts){
            std::cout<<"Partition #"<<pnumb++<<":";
            for (auto n:psg) std::cout<<" "<<(n? 1:0);
            std::cout<<std::endl;
        }
        std::cout<<std::endl;
        bool valid_part=true;
        if (parts.size()==2) {
            for (auto i = 0; i < bubbly.nodes.size(); ++i) {
                if (0 == i % 3) {
                    if (!parts[0][i] or !parts[1][i]) valid_part = false;
                } else if (parts[0][i] == parts[1][i]) valid_part = false;
                if (1 == i%3) if (parts[0][i]==parts[0][i+1]) valid_part = false;
            }
        }
        else valid_part=false;
        if (valid_part) {
            auto parts_score = partitioner.score_partition_set(bubbly, parts, tp);
            auto subgraphs = partitioner.partitions_as_subgraphs(bubbly, parts);
            std::cout << "Partition solutions as nodes:" << std::endl;
            pnumb = 0;
            for (auto &psg:subgraphs) {
                std::cout << "Partition #" << pnumb++ << ":";
                for (auto n:psg.nodes) std::cout << " " << n;
                std::cout << std::endl;

                sg.join_path(SequenceGraphPath(sg,psg.nodes));
            }
            std::cout << "Scores: " << parts_score.first << " " << parts_score.second << std::endl;
            ++solved_count;
        }
        else {
            std::cout<<"No valid solution found"<<std::endl;

        }
        //TODO: check the partition subgraphs can create paths, create paths if so, solve region.
        //TODO: find regions that are not "bubbly" but "repeaty".
    }
    std::cout<<"Bubbly paths solved trivially: "<<solved_count<<" / "<<bubblies.size()<<std::endl;
    //scaff.expand_bubbly_subgraphs();
    //scaff.pop_unsupported_shortbubbles();
    /*sg.join_all_unitigs();
    for (auto &m:mappers) {
        std::cout<<"removing obsolete mappings from "<<m.read1filename<<" and "<<m.read2filename<<std::endl;
        m.remove_obsolete_mappings();
        m.print_stats();
        m.remap_reads();
        m.print_stats();
    }*/
    //TODO: a lot of repeats are small repeats creating a big "loop", account for those!
    //scaff.find_canonical_repeats();






    sg.write_to_gfa(output_prefix+"_scaffolded.gfa");
    return 0;
}

