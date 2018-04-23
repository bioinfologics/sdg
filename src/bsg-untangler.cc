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
    bool repeat_expansion=false, neighbour_connection=false,neighbour_connection_graph=false, bubbly_paths=false;
    bool unroll_loops=false,pop_errors=false;
    bool explore_homopolymers=false;
    try
    {
        cxxopts::Options options("bsg-untangler", "graph-based repeat resolution and haplotype separation");

        options.add_options()
                ("help", "Print help")
                ("w,workspace", "input workspace", cxxopts::value<std::string>(workspace_file))
                ("o,output", "output file prefix", cxxopts::value<std::string>(output_prefix))
                ("l,unroll_loops", "unroll simple loops",cxxopts::value<bool>(unroll_loops))
                ("e,pop_errors", "pop unsupported short-bubbles (as errors)",cxxopts::value<bool>(pop_errors))
                ("r,repeat_expansion","run tag-based repeat expansion", cxxopts::value<bool>(repeat_expansion))
                ("b,bubbly_paths","run bubbly paths phasing", cxxopts::value<bool>(bubbly_paths))
                ("n,neighbour_connection","run tag-based repeat neighbour_connection", cxxopts::value<bool>(neighbour_connection))
                ("c,neighbour_connection_graph","create a tag-imbalance-based neighbour gfa", cxxopts::value<bool>(neighbour_connection_graph))
                ("explore_homopolymers","explore_homopolymers (experimental)", cxxopts::value<bool>(explore_homopolymers))
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

    if (unroll_loops){
        Untangler u(ws);
        u.unroll_simple_loops();
        /*ws.sg.join_all_unitigs();
        ws.kci.reindex_graph();
        for (auto &m:ws.linked_read_mappers) {
            m.remap_all_reads();
        }*/
    }
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
        ws.sg.write_to_gfa(output_prefix+"_repeat_solved.gfa");
        ws.kci.reindex_graph();
        for (auto &m:ws.linked_read_mappers) {
            m.remap_all_reads();
        }
    }
    if (bubbly_paths){
        Untangler u(ws);
        u.solve_bubbly_paths();
        ws.sg.join_all_unitigs();
        ws.sg.write_to_gfa(output_prefix+"_bubbly_solved.gfa");
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
    if (explore_homopolymers){
        //get all bubbles that are short
        Untangler u(ws);
        for (auto b:u.find_bubbles(360,401)){
            auto tA=ws.linked_read_mappers[0].get_node_tags(b.first);
            auto tB=ws.linked_read_mappers[0].get_node_tags(b.second);
            std::set<bsg10xTag> shared;
            std::set_intersection(tA.begin(),tA.end(),tB.begin(),tB.end(),std::inserter(shared,shared.end()));
            if (shared.empty()) continue;
            std::cout<<"Bubble analysed: " << b.first << " and " << b.second << " share "<< shared.size()<< " tags" <<std::endl;
            for (auto t: shared) std::cout<<"Tag "<<t<<" has "<<ws.linked_read_datastores[0].get_tag_reads(t).size()<<std::endl;
            std::set<uint64_t> readsA,readsB;
            std::cout<<"Reads in shared tags in A: ";
            for (auto rm:ws.linked_read_mappers[0].reads_in_node[llabs(b.first)]) {
                readsA.insert(rm.read_id);
                if (shared.count(ws.linked_read_datastores[0].get_read_tag(rm.read_id))>0) std::cout<<" "<<rm.read_id;
            }
            std::cout<<std::endl;
            std::cout<<"Reads in shared tags in B: ";
            for (auto rm:ws.linked_read_mappers[0].reads_in_node[llabs(b.second)]) {
                readsB.insert(rm.read_id);
                if (shared.count(ws.linked_read_datastores[0].get_read_tag(rm.read_id))>0) std::cout<<" "<<rm.read_id;
            }
            std::cout<<std::endl;

            for (auto ra:readsA){
                auto rb=(ra%2==1 ? ra+1 : ra-1);
                if (readsB.count(rb)>0) {
                    std::cout<<"Read "<<ra<<" in A and read"<<rb<<" in B!!!"<<std::endl;
                    std::cout<<ws.linked_read_datastores[0].get_read_sequence(ra)<<std::endl;
                    std::cout<<ws.linked_read_datastores[0].get_read_sequence(rb)<<std::endl;
                }
            }
        }
        //tags shared tags

        //reads from the shared tags
    }
    ws.dump_to_disk(output_prefix+"_untangled.bsgws");
    return 0;
}

