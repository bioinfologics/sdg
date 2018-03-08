#include <iostream>
#include <fstream>
#include <sglib/WorkSpace.hpp>
#include <sglib/processors/Untangler.hpp>
#include "sglib/logger/OutputLog.h"
#include "cxxopts.hpp"

int main(int argc, char * argv[]) {
    std::cout << "Welcome to bsg-untangler"<<std::endl<<std::endl;
    std::cout << "Git origin: " << GIT_ORIGIN_URL << " -> "  << GIT_BRANCH << std::endl;
    std::cout << "Git commit: " << GIT_COMMIT_HASH << std::endl<<std::endl;
    std::cout << "Executed command:"<<std::endl;
    for (auto i=0;i<argc;i++) std::cout<<argv[i]<<" ";
    std::cout<<std::endl<<std::endl;

    std::string workspace_file,output_prefix;
    bool stats_only=0,reindex_kci=0;
    uint64_t max_mem_gb=4;
    sglib::OutputLogLevel=sglib::LogLevels::DEBUG;
    bool repeats_first=1,pop_bubles=0,skip_remap=0,print_HSPNPs=1,haplotype_walk=0;
    std::string verbose_log="";
    try
    {
        cxxopts::Options options("bsg-untangler", "graph-based repeat resolution and haplotype separation");

        options.add_options()
                ("help", "Print help")
                ("w,workspace", "input workspace", cxxopts::value<std::string>(workspace_file))
                ("o,output", "output file prefix", cxxopts::value<std::string>(output_prefix))
                ("max_mem", "maximum_memory when mapping (GB, default: 4)", cxxopts::value<uint64_t>(max_mem_gb));

        options.add_options("Heuristics")
                ("skip_remap", "skip remapping reads on loading workspace (default false)", cxxopts::value<bool>(skip_remap))
                ("pop_short", "pop unsupported short bubbles first (default false)", cxxopts::value<bool>(pop_bubles))
                ("repeats_first", "solve_repeats before bubble phasing (default true)", cxxopts::value<bool>(repeats_first))
                ("print_hspnps", "print HSPNPs (default true)", cxxopts::value<bool>(print_HSPNPs))
                ("hwalk", "haplotype walk (default false)", cxxopts::value<bool>(haplotype_walk))
                ("heuristics_verbose_log", "dump heuristics verbose log to file", cxxopts::value<std::string>(verbose_log))
                ;



        auto result(options.parse(argc, argv));

        if (result.count("help"))
        {
            std::cout << options.help({"","Heuristics"}) << std::endl;
            exit(0);
        }

        if (result.count("w")!=1 or result.count("o")!=1) {
            throw cxxopts::OptionException(" please specify input workspace and output prefix");
        }



    } catch (const cxxopts::OptionException& e)
    {
        std::cout << "Error parsing options: " << e.what() << std::endl << std::endl
                <<"Use option --help to check command line arguments." << std::endl;
        exit(1);
    }

    std::cout<<std::endl;
    WorkSpace ws;
    sglib::OutputLog()<<"Loading Workspace..."<<std::endl;
    ws.load_from_disk(workspace_file);
    ws.add_log_entry("bsg-untangler run started");
    sglib::OutputLog()<<"Loading Workspace DONE"<<std::endl;
    ws.kci.compute_compression_stats();
    for (auto &m:ws.linked_read_mappers) m.memlimit=max_mem_gb*1024*1024*1024;
    if (skip_remap==false) {
        sglib::OutputLog()<<"Mapping reads..."<<std::endl;
        for (auto &m:ws.linked_read_mappers) {
            m.update_graph_index();
            m.map_reads();
        }
        ws.add_log_entry("reads re-mapped to current graph");
        ws.dump_to_disk(output_prefix+"_mapped.bsgws");
        sglib::OutputLog()<<"Mapping reads DONE."<<std::endl;
    }




    std::ofstream verbose_log_file;
    if (verbose_log!=""){
        verbose_log_file.open(verbose_log);
    }

    if (repeats_first) {
        std::cout << std::endl << "Solving trivial repeats" << std::endl;
        bool mod = true;
        int srpass = 0;
        while (mod) {
            Untangler u(ws);
            std::unordered_set<uint64_t> reads_to_remap;
            mod = false;
            ++srpass;
            uint64_t solved_paths=u.solve_canonical_repeats_by_tags(reads_to_remap);
            if (solved_paths>0) {
                mod=true;
                ws.kci.reindex_graph();
                //ws.sg.write_to_gfa(output_prefix + "_solved_repeats_" + std::to_string(srpass) + ".gfa",unsolved_repeats);
                std::cout << reads_to_remap.size() << " pairs of read to remap" << std::endl;
                ws.linked_read_mappers[0].update_graph_index();
                ws.linked_read_mappers[0].map_reads(reads_to_remap);
                ws.add_log_entry(std::to_string(solved_paths)+" new paths by solving canonical repeats, "+
                                         std::to_string(reads_to_remap.size()) +"reads re-mapped");
            }

        }
    }

    Untangler u(ws);
    if (haplotype_walk){
        u.extend_HSPNPs_by_tagwalking();
    }
    if (print_HSPNPs) {
        std::cout<<std::endl<<"Finding HSPNPs"<<std::endl;
        auto hps=u.get_all_HSPNPs();
        std::cout << "Starting with " << hps.size() << " HSPNPs" << std::endl;
        std::ofstream hpsfile(output_prefix + "_hps.csv");
        auto hspnp_id = 1;
        for (auto hp:hps) { //TODO: print node numbers and alternative tags -> then cluster/find neighbours
            std::vector<prm10xTag_t> h1tags, h2tags;
            for (auto rm:ws.linked_read_mappers[0].reads_in_node[(hp.first > 0 ? hp.first : -hp.first)])
                h1tags.push_back(ws.linked_read_mappers[0].datastore.get_read_tag(rm.read_id));
            for (auto rm:ws.linked_read_mappers[0].reads_in_node[(hp.second > 0 ? hp.second : -hp.second)])
                h2tags.push_back(ws.linked_read_mappers[0].datastore.get_read_tag(rm.read_id));
            std::sort(h1tags.begin(), h1tags.end());
            std::sort(h2tags.begin(), h2tags.end());
            hpsfile << hspnp_id << ",A," << hp.first;
            for (auto t:h1tags) if (t) hpsfile << "," << t;
            hpsfile << std::endl;
            hpsfile << hspnp_id << ",B," << hp.second;
            for (auto t:h2tags) if (t) hpsfile << "," << t;
            hpsfile << std::endl;
            ++hspnp_id;
        }
    }


    sglib::OutputLog()<<"Dumping final Workspace..."<<std::endl;
    ws.dump_to_disk(output_prefix+"_final.bsgws");
    sglib::OutputLog()<<"Dumping scaffolded GFA..."<<std::endl;
    ws.sg.write_to_gfa(output_prefix+"_scaffolded.gfa");
    sglib::OutputLog()<<"All DONE!!!"<<std::endl;
    return 0;
}

