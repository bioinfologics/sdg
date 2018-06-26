#include <iostream>
#include <fstream>
#include <sglib/WorkSpace.hpp>
#include <sglib/processors/Untangler.hpp>
#include <sglib/processors/FlowFollower.hpp>
#include <sglib/processors/LinkageUntangler.hpp>
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
    bool repeat_expansion=false, neighbour_connection_graph=false, bubbly_paths=false;
    bool unroll_loops=false,pop_errors=false, paired_scaff=false, select_hspnps=false;
    bool explore_homopolymers=false;
    uint64_t min_backbone_node_size=1000;
    float min_backbone_ci=0.5;
    float max_backbone_ci=1.25;
    float tag_imbalance_ends=.1;
    int min_pairs=7;
    int min_shared_tags=10;
    try
    {
        cxxopts::Options options("bsg-untangler", "graph-based repeat resolution and haplotype separation");

        options.add_options()
                ("help", "Print help")
                ("w,workspace", "input workspace", cxxopts::value<std::string>(workspace_file))
                ("o,output", "output file prefix", cxxopts::value<std::string>(output_prefix))
                ("p,paired_reads","paired read scaffolding (experimental)",cxxopts::value<bool>(paired_scaff))
                ("l,unroll_loops", "unroll simple loops",cxxopts::value<bool>(unroll_loops))
                ("e,pop_errors", "pop unsupported short-bubbles (as errors)",cxxopts::value<bool>(pop_errors))
                ("r,repeat_expansion","run tag-based repeat expansion", cxxopts::value<bool>(repeat_expansion))
                ("b,bubbly_paths","run bubbly paths phasing", cxxopts::value<bool>(bubbly_paths))
                //("n,neighbour_connections","run tag-based repeat neighbour_connection", cxxopts::value<bool>(neighbour_connection))
                ("c,neighbour_connection_backbones","create tag-imbalance-based neighbour backbones", cxxopts::value<bool>(neighbour_connection_graph))
                ("select_hspnp","select only Haplotype Specific Parallel Node Pairs to base scaffolding", cxxopts::value<bool>(select_hspnps))
                ("min_backbone_node_size","minimum size of nodes to use in backbones",cxxopts::value<uint64_t>(min_backbone_node_size))
                ("min_pairs","minimum number of pairs to connect two nodes",cxxopts::value<int>(min_pairs))
                ("min_backbone_ci","minimum ci of nodes to use in backbones",cxxopts::value<float>(min_backbone_ci))
                ("max_backbone_ci","minimum ci of nodes to use in backbones",cxxopts::value<float>(max_backbone_ci))
                ("tag_imbalance_ends","percentage of node to use as tag-imbalanced end",cxxopts::value<float>(tag_imbalance_ends))
                ("min_shared_tags","minimum shared tags to evaluate tag-imbalanced connection",cxxopts::value<int>(min_shared_tags))
                //("explore_homopolymers","explore_homopolymers (experimental)", cxxopts::value<bool>(explore_homopolymers))
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
    if (paired_scaff){
        if (!ws.linked_read_mappers.empty()) {
            for (auto round=1;round<11;++round) {
                sglib::OutputLog()<<"STARTING ROUND #"<<std::to_string(round)<<std::endl;
                LinkageUntangler lu(ws);
                if (select_hspnps) lu.select_nodes_by_HSPNPs(min_backbone_node_size,min_backbone_ci,max_backbone_ci);
                else lu.select_nodes_by_size_and_ci(min_backbone_node_size,min_backbone_ci,max_backbone_ci);

                std::unordered_set<sgNodeID_t> selnodes;
                for (sgNodeID_t n=1;n<ws.sg.nodes.size();++n) if (lu.selected_nodes[n]) selnodes.insert(n);
                lu.report_node_selection();
                /*auto topology_ldg=lu.make_topology_linkage(10);
                topology_ldg.report_connectivity();
                ws.sg.write_to_gfa("topology_links.gfa",{},{},{},topology_ldg.links);*/
                LinkageDiGraph gldg(ws.sg);
                /*if (!ws.paired_read_mappers.empty()) {
                    auto pair_ldg = lu.make_paired_linkage(min_pairs);
                    pair_ldg.report_connectivity();
                    gldg.add_links(pair_ldg);
                    ws.sg.write_to_gfa("pair_links.gfa", {}, {}, {}, pair_ldg.links);
                    pair_ldg.remove_transitive_links(10);
                    pair_ldg.report_connectivity();
                    ws.sg.write_to_gfa("pair_links_no_transitive.gfa", {}, {}, {}, pair_ldg.links);
                }*/
                sglib::OutputLog()<<"Eliminating N-N nodes..."<<std::endl;
                //HACK: eliminate nodes with N-N and try again.
                auto pre_tag_ldg = lu.make_tag_linkage(min_shared_tags);
                pre_tag_ldg.remove_transitive_links(10);
                pre_tag_ldg.report_connectivity();
                uint64_t remNN=0;
                for (auto n=1;n<ws.sg.nodes.size();++n){
                    if (lu.selected_nodes[n]){
                        if (pre_tag_ldg.get_fw_links(n).size()>1 and pre_tag_ldg.get_bw_links(n).size()>1) {
                            lu.selected_nodes[n]=false;
                            ++remNN;
                        }
                    }
                }
                sglib::OutputLog()<<"Re-trying tag connection after eliminating "<<remNN<<" N-N nodes"<<std::endl;
                auto tag_ldg = lu.make_tag_linkage(min_shared_tags);
                tag_ldg.remove_transitive_links(10);
                tag_ldg.report_connectivity();
                ws.sg.write_to_gfa(output_prefix + "_tag_nt_" + std::to_string(round) + ".gfa", {}, {}, selnodes, tag_ldg.links);
                sglib::OutputLog() << "Simplifying linear paths" << std::endl;
                lu.expand_linear_regions(tag_ldg);
                auto joined=ws.sg.join_all_unitigs();
                ws.sg.write_to_gfa(output_prefix + "_after_expansion_" + std::to_string(round) + ".gfa", {}, {});
                //sglib::OutputLog()<<"TODO: remap reads and re-start the whole thing..."<<std::endl;
                ws.kci.reindex_graph();
                ws.remap_all();
                if (joined==0) break;
            }
            /*
            sglib::OutputLog()<<"Eliminating N-N nodes..."<<std::endl;
            //HACK: eliminate nodes with N-N and try again.
            auto sel_orig=lu.selected_nodes;
            uint64_t remNN=0;
            for (auto n=1;n<ws.sg.nodes.size();++n){
                if (lu.selected_nodes[n]){
                    if (tag_ldg.get_fw_links(n).size()>1 and tag_ldg.get_bw_links(n).size()>1) {
                        lu.selected_nodes[n]=false;
                        ++remNN;
                    }
                }
            }
            sglib::OutputLog()<<"Re-trying tag connection after eliminating "<<remNN<<" N-N nodes"<<std::endl;
            auto tag_ldg_noNN = lu.make_tag_linkage(min_shared_tags);
            tag_ldg_noNN.report_connectivity();
            tag_ldg_noNN.remove_transitive_links(10);
            tag_ldg_noNN.report_connectivity();
            ws.sg.write_to_gfa(output_prefix+"_tag_noNN_nt.gfa", {}, {}, selnodes, tag_ldg_noNN.links);
            lu.selected_nodes=sel_orig;*/



            /*exit(0);

            tag_ldg.report_connectivity();
            gldg.add_links(tag_ldg);
            ws.sg.write_to_gfa("tag_links.gfa", {}, {}, {}, tag_ldg.links);
            auto tag_hspnp_ldg=lu.filter_linkage_to_hspnp_duos(min_backbone_node_size,min_backbone_ci,max_backbone_ci,tag_ldg);
            tag_hspnp_ldg.report_connectivity();
            //ws.sg.write_to_gfa("tag_links_hspnps_coherent_selonly.gfa", {}, {},selnodes, tag_hspnp_ldg.links);
            tag_hspnp_ldg.remove_transitive_links(10);
            tag_hspnp_ldg.report_connectivity();
            //ws.sg.write_to_gfa("tag_links_hspnps_coherent_selonly_no_transitive.gfa", {}, {}, selnodes, tag_hspnp_ldg.links);
            lu.expand_trivial_repeats(tag_hspnp_ldg);
            ws.sg.join_all_unitigs();
            ws.sg.write_to_gfa("after_trivial_hspnp_expansion.gfa");//, {}, {},{}, tag_hspnp_ldg.links);
            //TODO: remap reads!
            ws.sg.create_index();
            for (auto &m:ws.paired_read_mappers) {
                sglib::OutputLog()<<"Mapping reads from paired library..."<<std::endl;
                m.map_reads();
                m.print_stats();
                ws.add_log_entry("reads from "+m.datastore.filename+" re-mapped to current graph");
                sglib::OutputLog()<<"Mapping reads from paired library DONE."<<std::endl;
            }
            for (auto &m:ws.linked_read_mappers) {
                sglib::OutputLog()<<"Mapping reads from linked library..."<<std::endl;
                m.map_reads();
                ws.add_log_entry("reads from "+m.datastore.filename+" re-mapped to current graph");
                sglib::OutputLog()<<"Mapping reads from linked library DONE."<<std::endl;
            }

            tag_ldg.remove_transitive_links(10);
            tag_ldg.report_connectivity();
            ws.sg.write_to_gfa("tag_links_no_transitive.gfa", {}, {}, {}, tag_ldg.links);
            ws.sg.write_to_gfa("tag_links_no_transitive_selected_only.gfa", {}, {}, selnodes, tag_ldg.links);
            //HACK: eliminate nodes with N-N and try again.
            auto sel_orig=lu.selected_nodes;
            uint64_t remNN=0;
            for (auto n=1;n<ws.sg.nodes.size();++n){
                if (lu.selected_nodes[n]){
                    if (tag_ldg.get_fw_links(n).size()>1 and tag_ldg.get_bw_links(n).size()>1) {
                        lu.selected_nodes[n]=false;
                        ++remNN;
                    }
                }
            }
            sglib::OutputLog()<<"Re-trying tag connection after eliminating "<<remNN<<" N-N nodes"<<std::endl;
            auto tag_ldg_noNN = lu.make_tag_linkage(min_shared_tags);
            tag_ldg_noNN.report_connectivity();
            //gldg.add_links(tag_ldg);
            ws.sg.write_to_gfa("tag_links_noNN.gfa", {}, {}, {}, tag_ldg_noNN.links);
            tag_ldg_noNN.remove_transitive_links(10);
            tag_ldg_noNN.report_connectivity();
            ws.sg.write_to_gfa("tag_links_noNN_no_transitive.gfa", {}, {}, {}, tag_ldg_noNN.links);
            ws.sg.write_to_gfa("tag_links_noNN_no_transitive_selected_only.gfa", {}, {}, selnodes, tag_ldg_noNN.links);
            lu.selected_nodes=sel_orig;*/


        }
        /*gldg.report_connectivity();
        gldg.add_links(gldg);
        ws.sg.write_to_gfa("general_links.gfa", {}, {}, {}, gldg.links);
        gldg.remove_transitive_links(10);
        gldg.report_connectivity();
        ws.sg.write_to_gfa("general_links_no_transitive.gfa", {}, {}, {}, gldg.links);*/
        //PairedReadLinker prl(ws,u);
        //prl.generate_links_size_ci(min_backbone_node_size,min_backbone_ci,max_backbone_ci,5);
        //prl.generate_links_hspnp();
        //ws.sg.write_to_gfa("prl_hspnp_links.gfa",{},{},{},prl.links);
        //std::cout<<"calling remove_transitive_links"<<std::endl;
        //prl.remove_transitive_links();
        //std::cout<<"remove_transitive_links finished"<<std::endl;
        /*for (auto lp:prl.find_local_problems(15000)){
            std::cout<<"Local problem, frontiers: ";
            for (auto n:lp) if (ws.sg.nodes[llabs(n)].sequence.size()>=15000) std::cout<<" "<<n;
            std::cout<<std::endl<<"                internal: ";
            for (auto n:lp) if (ws.sg.nodes[llabs(n)].sequence.size()<15000) std::cout<<" "<<n;
            std::cout<<std::endl<<"    ";
            for (auto n:lp) if (ws.sg.nodes[llabs(n)].sequence.size()) std::cout<<" seq"<<llabs(n)<<",";
            std::cout<<std::endl<<std::endl;
            prl.solve_local_problem(lp);
        }*/
    }
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
    if (neighbour_connection_graph){
        Untangler u(ws);
        auto backbones=u.create_backbones(min_backbone_node_size,min_backbone_ci, max_backbone_ci,tag_imbalance_ends,min_shared_tags);
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

