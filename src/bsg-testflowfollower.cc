#include <iostream>
#include <fstream>
#include <sglib/WorkSpace.hpp>
#include <sglib/processors/Untangler.hpp>
#include <sglib/processors/FlowFollower.hpp>
#include "sglib/logger/OutputLog.h"
#include "cxxopts.hpp"

int main(int argc, char * argv[]) {
    std::cout << "Welcome to bsg-testflowfollower"<<std::endl<<std::endl;
    std::cout << "Git origin: " << GIT_ORIGIN_URL << " -> "  << GIT_BRANCH << std::endl;
    std::cout << "Git commit: " << GIT_COMMIT_HASH << std::endl<<std::endl;
    std::cout << "Executed command:"<<std::endl;
    for (auto i=0;i<argc;i++) std::cout<<argv[i]<<" ";
    std::cout<<std::endl<<std::endl;

    std::string workspace_file,output_prefix;
    uint64_t max_mem_gb=4;
    sglib::OutputLogLevel=sglib::LogLevels::DEBUG;
    sgNodeID_t start_node=0,end_node=0;
    uint64_t max_nodes=10000;
    uint16_t select_min_tags=0,select_max_tags=0;
    float select_min_ci=.75,select_max_ci=1.5;
    uint32_t select_min_size=399,select_max_size=5000;
    try
    {
        cxxopts::Options options("bsg-testflowfollower", "a test following flows in a parallel region of the graph");

        options.add_options()
                ("help", "Print help")
                ("w,workspace", "input workspace", cxxopts::value<std::string>(workspace_file))
                ("o,output", "output file prefix", cxxopts::value<std::string>(output_prefix))
                //("s,start", "starting node (with direction flowing towards end)", cxxopts::value<sgNodeID_t>(start_node))
                //("e,end", "ending node (with direction flowing from start)", cxxopts::value<sgNodeID_t>(end_node))
                ("s,min_size", "node selection, min size", cxxopts::value<uint32_t>(select_min_size))
                ("S,max_size", "node selection, max size", cxxopts::value<uint32_t>(select_max_size))
                ("c,min_ci", "node selection, min compression index", cxxopts::value<float>(select_min_ci))
                ("C,max_ci", "node selection, max compression index", cxxopts::value<float>(select_max_ci))
                ("t,min_tags", "node selection, min tags", cxxopts::value<uint16_t>(select_min_tags))
                ("T,max_tags", "node selection, max tags", cxxopts::value<uint16_t>(select_max_tags));
                //("n,max_nodes", "maximum number of nodes to explore (default: 10000)", cxxopts::value<uint64_t >(max_nodes));



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
    //ws.sg.write_to_gfa("initial_graph.gfa",{},{},{});
    ws.linked_read_datastores[0].dump_tag_occupancy_histogram("tag_occupancy.csv");
    ws.add_log_entry("bsg-testflowfollower run started");
    sglib::OutputLog()<<"Loading Workspace DONE"<<std::endl;
    //ws.kci.reindex_graph();
    ws.kci.compute_compression_stats();
    //for (auto &m:ws.linked_read_mappers) m.memlimit=max_mem_gb*1024*1024*1024;
//TODO: move region selection into SG
//    std::unordered_set<sgNodeID_t> region_nodes;
//    if (select_min_tags==0) {
//        if (start_node != 0 and end_node != 0) {
//            //find the path
//            region_nodes.insert(start_node);
//            uint64_t last_size = 0;
//            bool end_found = false;
//            while (region_nodes.size() < max_nodes and last_size < region_nodes.size()) {
//                last_size = region_nodes.size();
//                std::vector<sgNodeID_t> new_nodes;
//                for (auto n:region_nodes) {
//                    if (n != end_node) {
//                        for (auto l:ws.sg.get_fw_links(n)) new_nodes.push_back(l.dest);
//                    } else end_found = true;
//                    if (n != start_node) {
//                        for (auto l:ws.sg.get_bw_links(n)) new_nodes.push_back(-l.dest);;
//                    }
//                }
//                for (auto n:new_nodes) region_nodes.insert(n);
//            }
//            if (region_nodes.size() >= max_nodes) {
//                std::cout << "Aborting after growing the region to " << region_nodes.size() << " nodes" << std::endl;
//                exit(0);
//            }
//            std::cout << "Found a region with " << region_nodes.size() << " nodes" << std::endl;
//
//            ws.sg.write_to_gfa(output_prefix + "_flowregion.gfa", {}, {}, region_nodes);
//        } else {
//            ws.kci.reindex_graph();
//            ws.kci.compute_compression_stats();
//            Untangler untangler(ws);
//            auto HSPNPs = untangler.get_all_HSPNPs();
//            for (auto h:HSPNPs) {
//                region_nodes.insert(h.first);
//                region_nodes.insert(h.second);
//            }
//            std::cout << "Using the whole graph with " << region_nodes.size() << " nodes from HSPNPs" << std::endl;
//        }
//    }

    FlowFollower ff(ws,{});
    //TODO: add option to include HSPNPs
//    ws.kci.reindex_graph();
//            ws.kci.compute_compression_stats();
//            Untangler untangler(ws);
//            auto HSPNPs = untangler.get_all_HSPNPs();
//            for (auto h:HSPNPs) {
//                region_nodes.insert(h.first);
//                region_nodes.insert(h.second);
        ws.add_log_entry("Creating paths from linked reads, starting nodes: "+
                         std::to_string(select_min_size)+ "-" + std::to_string(select_max_size) + " bp "+
                         std::to_string(select_min_tags)+"-"+std::to_string(select_min_tags)+ " tags " +
                         std::to_string(select_min_ci)+ "-" + std::to_string(select_max_ci) + " CI");
        ff.select_from_all_nodes(select_min_size,select_max_size,select_min_tags, select_max_tags,
                                 select_min_ci,select_max_ci);
        //ff.create_flows();
        ff.create_flows_fast();

//    }
    sglib::OutputLog()<<"Dumping final Workspace..."<<std::endl;
    ws.dump_to_disk(output_prefix+"_final.bsgws");
//    sglib::OutputLog()<<"Dumping scaffolded GFA..."<<std::endl;
//    ws.sg.write_to_gfa(output_prefix+"_scaffolded.gfa");
//    sglib::OutputLog()<<"All DONE!!!"<<std::endl;
//    return 0;
}

