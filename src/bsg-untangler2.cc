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
    bool repeat_expansion=false, neighbour_connection=false;
    try
    {
        cxxopts::Options options("bsg-untangler", "graph-based repeat resolution and haplotype separation");

        options.add_options()
                ("help", "Print help")
                ("w,workspace", "input workspace", cxxopts::value<std::string>(workspace_file))
                ("o,output", "output file prefix", cxxopts::value<std::string>(output_prefix))
                ("r,repeat_expansion","run tag-based repeat expansion", cxxopts::value<bool>(repeat_expansion))
                ("n,neighbour_connection","run tag-based repeat neighbour_connection", cxxopts::value<bool>(neighbour_connection))
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
    if (ws.path_datastores.size()==0) {
        sglib::OutputLog()<<"Finishing early because there's no path_datastores[0]"<<std::endl;
        return 1;
    }

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
    if (repeat_expansion) {
        //==================== Development code (i.e. random tests!) ==============
        Untangler u(ws);
        u.expand_canonical_repeats_by_tags(.5,1.5);
        if (!ws.sg.is_sane()) {
            sglib::OutputLog()<<"ERROR: sg.is_sane() = false"<<std::endl;
            return 1;
        }
        ws.sg.join_all_unitigs();
        ws.kci.reindex_graph();
        for (auto &m:ws.linked_read_mappers) {
            m.remap_all_reads();
        }
        ws.dump_to_disk(output_prefix+"_repeats_expanded.bsgws");


    }
    if (neighbour_connection){
        Untangler u(ws);
        u.connect_neighbours();
    }
    return 0;
}

