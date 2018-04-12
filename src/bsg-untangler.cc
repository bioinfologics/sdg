#include <iostream>
#include <fstream>
#include <sglib/WorkSpace.hpp>
#include <sglib/processors/Untangler.hpp>
#include <sglib/processors/FlowFollower.hpp>
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
    sglib::OutputLogLevel=sglib::LogLevels::DEBUG;
    bool repeat_expansion=false, neighbour_connection=false,neighbour_connection_graph=false, bubbly_paths=false, pop_errors=false;
    try
    {
        cxxopts::Options options("bsg-untangler", "graph-based repeat resolution and haplotype separation");

        options.add_options()
                ("help", "Print help")
                ("w,workspace", "input workspace", cxxopts::value<std::string>(workspace_file))
                ("o,output", "output file prefix", cxxopts::value<std::string>(output_prefix))
                ("e,pop_errors", "pop unsupported short-bubbles (as errors)",cxxopts::value<bool>(pop_errors))
                ("r,repeat_expansion","run tag-based repeat expansion", cxxopts::value<bool>(repeat_expansion))
                ("b,bubbly_paths","run bubbly paths phasing", cxxopts::value<bool>(bubbly_paths))
                ("n,neighbour_connection","run tag-based repeat neighbour_connection", cxxopts::value<bool>(neighbour_connection))
                ("c,neighbour_connection_graph","create a tag-imbalance-based neighbour gfa", cxxopts::value<bool>(neighbour_connection_graph))
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

    /*if (!devel_code) {

        Untangler u(ws);

        u.analise_paths_through_nodes();
        FlowFollower ff(ws);

        ff.load_flows_from_paths(ws.path_datastores.back());
        std::ofstream skf("skated_paths.fasta");
        for (auto sp:ff.skate_from_all(4, 10000)) {
            skf << ">" << sp.nodes[0] << "_" << sp.nodes.back() << std::endl << sp.get_sequence() << std::endl;
        }

        std::cout << "TODO: pop error-bp bubbles by kci and paths" << std::endl;
        u.pop_errors_by_ci_and_paths();
        std::cout
                << "TODO: Solve bubbly paths by kci (both local and bubly-path-total) and paths through collapses (expand without joins, update paths!)"
                << std::endl;
        std::cout << "TODO: Join unitigs and remap all reads" << std::endl;
        std::cout << "TODO: skate through crap, generate a copy of the solution and re-connect" << std::endl;
        std::cout << "TODO: Join unitigs and remap all reads" << std::endl;
        sglib::OutputLog() << "All DONE!!!" << std::endl;

        //u.analise_paths_through_nodes();
    }*/
    if (pop_errors){
        Untangler u(ws);
        u.pop_errors_by_ci_and_paths();
        ws.sg.join_all_unitigs();
        ws.kci.reindex_graph();
        for (auto &m:ws.linked_read_mappers) {
            m.remap_all_reads();
        }
    }
    if (repeat_expansion) {
        //==================== Development code (i.e. random tests!) ==============
        Untangler u(ws);
        u.expand_canonical_repeats_by_tags(.5,1.25,20);
        if (!ws.sg.is_sane()) {
            sglib::OutputLog()<<"ERROR: sg.is_sane() = false"<<std::endl;
            return 1;
        }
        ws.sg.join_all_unitigs();
        ws.kci.reindex_graph();
        for (auto &m:ws.linked_read_mappers) {
            m.remap_all_reads();
        }
    }
    if (bubbly_paths){
        Untangler u(ws);
        u.solve_bubbly_paths();
        ws.sg.join_all_unitigs();
        //ws.sg.write_to_gfa(output_prefix+"_bubbly_solved.gfa");
        ws.kci.reindex_graph();
        for (auto &m:ws.linked_read_mappers) {
            m.remap_all_reads();
        }

    }
    if (neighbour_connection){
        Untangler u(ws);
        int i=0;
        for (auto minsize:{500,1000,1500,2500,5000,7500,10000}) {
            uint64_t last=1;
            while (last) {
                last = u.connect_neighbours(minsize, .5, 1.25, 50000);
                if (last) {
                    ws.sg.join_all_unitigs();
                    ws.kci.reindex_graph();
                    ws.sg.write_to_gfa(output_prefix + "_partial" + std::to_string(++i) + ".gfa");
                    for (auto &m:ws.linked_read_mappers) {
                        m.remap_all_reads();
                    }
                    //ws.dump_to_disk(output_prefix+"_partial"+std::to_string(++i)+".bsgws");
                }
            }
        }
        /*ws.kci.reindex_graph();
        for (auto &m:ws.linked_read_mappers) {
            m.remap_all_reads();
        }*/
    }
    if (neighbour_connection_graph){
        Untangler u(ws);
        auto backbones=u.create_backbones(.5, 1.25,.25);
        //auto tni=u.find_tag_neighbours_with_imbalance(5000, .5, 1.25,.25);
        //create a gfa with all nodes in nti, and connect them, dump the gfa.

        /*ws.kci.reindex_graph();
        for (auto &m:ws.linked_read_mappers) {
            m.remap_all_reads();
        }*/
    }
    ws.dump_to_disk(output_prefix+"_untangled.bsgws");
    return 0;
}

