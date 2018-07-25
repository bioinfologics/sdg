#include <iostream>
#include <fstream>
#include <sglib/WorkSpace.hpp>
#include <sglib/processors/Untangler.hpp>
#include <sglib/processors/FlowFollower.hpp>
#include <sglib/processors/LinkageUntangler.hpp>
#include <sglib/processors/LocalHaplotypeAssembler.hpp>
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
    int dev_max_lines=0;
    uint8_t dev_local_k=63;
    int dev_local_min_cvg=7;
    uint64_t dev_min_nodes=2,dev_min_total_size=0;
    std::string dev_create_linkage,dev_skate_linkage,dev_local_assembly_linkage;
    bool dev_local_patching=false;
    bool remap_reads=true;
    bool dump_gfa=false;
    bool dev_linkage_paths=false;
    std::string dev_linkage_stats;
    try
    {
        cxxopts::Options options("bsg-untangler", "graph-based haplotype separation");

        options.add_options()
                ("help", "Print help")
                ("w,workspace", "input workspace", cxxopts::value<std::string>(workspace_file))
                ("o,output", "output file prefix", cxxopts::value<std::string>(output_prefix))
                ("m,remap_reads","remap all reads on final workspace (default: true)",cxxopts::value<bool>(remap_reads))
                ("dump_gfa","dump the final graph in GFA (default: false)",cxxopts::value<bool>(dump_gfa));

        options.add_options("Backbone (unique anchors linkage)")
                ("min_backbone_node_size","minimum size of nodes to use in backbones",cxxopts::value<uint64_t>(min_backbone_node_size))
                ("min_backbone_ci","minimum ci of nodes to use in backbones",cxxopts::value<float>(min_backbone_ci))
                ("max_backbone_ci","minimum ci of nodes to use in backbones",cxxopts::value<float>(max_backbone_ci))
                ("select_hspnp","select only Haplotype Specific Parallel Node Pairs to base scaffolding", cxxopts::value<bool>(select_hspnps));

        options.add_options("Untangling functions")
                ("p,paired_reads","paired read scaffolding (experimental)",cxxopts::value<bool>(paired_scaff))
                ("l,unroll_loops", "unroll simple loops",cxxopts::value<bool>(unroll_loops))
                ("e,pop_errors", "pop unsupported short-bubbles (as errors)",cxxopts::value<bool>(pop_errors))
                ("r,repeat_expansion","run tag-based repeat expansion", cxxopts::value<bool>(repeat_expansion))
                ("b,bubbly_paths","run bubbly paths phasing", cxxopts::value<bool>(bubbly_paths))
                ("min_pairs","minimum number of pairs to connect two nodes",cxxopts::value<int>(min_pairs))
                ("tag_imbalance_ends","percentage of node to use as tag-imbalanced end",cxxopts::value<float>(tag_imbalance_ends))
                ("min_shared_tags","minimum shared tags to evaluate tag-imbalanced connection",cxxopts::value<int>(min_shared_tags));

        options.add_options("Development")
                ("dev_create_linkage","Creates and simplifies linkage and dumps to file",cxxopts::value<std::string>(dev_create_linkage))
                ("dev_linkage_paths", "tag linkage uses read pathing rather than simple mapping",cxxopts::value<bool>(dev_linkage_paths))
                ("dev_skate_linkage","Loads linkage from file and skates",cxxopts::value<std::string>(dev_skate_linkage))
                ("dev_local_assembly_linkage","Loads linkage from file and creates local assemblies",cxxopts::value<std::string>(dev_local_assembly_linkage))
                ("dev_linkage_stats","Loads linkage from file and computes local assemblies stats",cxxopts::value<std::string>(dev_linkage_stats))
                ("dev_max_lines","Limits lines to be skated on dev",cxxopts::value<int>(dev_max_lines))
                ("dev_min_nodes","Limits lines to be locally assembled on dev to at least min_nodes",cxxopts::value<uint64_t>(dev_min_nodes))
                ("dev_min_total_size","Limits lines to be locally assembled on dev to at least min_total_size",cxxopts::value<uint64_t>(dev_min_total_size))
                ("dev_local_k","k value for local assembly",cxxopts::value<uint8_t>(dev_local_k))
                ("dev_local_min_cvg","minimum coverga for local assemblies",cxxopts::value<int>(dev_local_min_cvg))
                ("dev_local_patching","run a whole round of linked lines and local patching, with a final remap", cxxopts::value<bool>(dev_local_patching));




        auto result(options.parse(argc, argv));

        if (result.count("help"))
        {
            std::cout << options.help({"","Backbone (unique anchors linkage)","Untangling functions","Development"}) << std::endl;
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

    //======= WORKSPACE LOAD AND CHECKS ======
    WorkSpace ws;
    sglib::OutputLog()<<"Loading Workspace..."<<std::endl;
    ws.load_from_disk(workspace_file);
    if (!ws.sg.is_sane()) {
        sglib::OutputLog()<<"ERROR: sg.is_sane() = false"<<std::endl;
        return 1;
    }
    //TODO: other checks? reads mapped to valid nodes and such?

    ws.add_log_entry("bsg-untangler run started");
    sglib::OutputLog()<<"Loading Workspace DONE"<<std::endl;

    //======= FUNCTIONS FOR REAL-LIFE USE ======

    if (unroll_loops){
        Untangler u(ws);
        u.unroll_simple_loops();
        ws.sg.join_all_unitigs();
    }
    if (pop_errors){
        Untangler u(ws);
        u.pop_errors_by_ci_and_paths(60,200);
        ws.sg.join_all_unitigs();
        ws.kci.reindex_graph();
        for (auto &m:ws.linked_read_mappers) {
            m.remap_all_reads();
        }
    }
    if (repeat_expansion) {
        //TODO: parameters are hardcoded here
        Untangler u(ws);
        u.expand_canonical_repeats_by_tags(.5,1.25,20);
        ws.sg.join_all_unitigs();
    }
    if (bubbly_paths){
        Untangler u(ws);
        u.solve_bubbly_paths();
        ws.sg.join_all_unitigs();
    }
    if (paired_scaff){
        sglib::OutputLog()<<"Creating node linkage from kmers..."<<std::endl;
        LinkageUntangler lu(ws);
        lu.make_paired_linkage_by_kmer(2,{0},false,true);
        exit(0);
        Untangler u(ws);
        sglib::OutputLog()<<"Popping errors..."<<std::endl;
        u.pop_errors_by_ci_and_paths(60,200);
        //auto juc=ws.sg.join_all_unitigs();
        //sglib::OutputLog()<<"Unitigs joined: "<<juc<<std::endl;
        //TODO: simple paired end repeat resolution, graph overlap expansion (alla arre), etc.

        lu.select_nodes_by_size_and_ci(1000,0,50);
        lu.report_node_selection();
        auto pld=lu.make_paired_linkage_pe(5);
        ws.sg.write_to_gfa(output_prefix+"_linkage.gfa", {}, {}, {}, pld.links);
        lu.expand_trivial_repeats(pld);
        auto juc=ws.sg.join_all_unitigs();
        sglib::OutputLog()<<"Unitigs joined: "<<juc<<std::endl;

    }

    //======= DEVELOPMENT FUNCTIONS ======
    if (!dev_create_linkage.empty()) {
        LinkageUntangler lu(ws);
        lu.select_nodes_by_size_and_ci(min_backbone_node_size,min_backbone_ci,max_backbone_ci);
        std::unordered_set<sgNodeID_t> selnodes;
        for (sgNodeID_t n=1;n<ws.sg.nodes.size();++n) if (lu.selected_nodes[n]) selnodes.insert(n);
        lu.report_node_selection();
        auto pre_tag_ldg = lu.make_tag_linkage(min_shared_tags,dev_linkage_paths);
        pre_tag_ldg.remove_transitive_links(10);
        pre_tag_ldg.report_connectivity();
        sglib::OutputLog()<<"Eliminating N-N nodes..."<<std::endl;
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
        auto tag_ldg = lu.make_tag_linkage(min_shared_tags,dev_linkage_paths);
        //maybe keep the full linkage to help the evaluation on simplification?
        auto full_tag_ldg=tag_ldg;
        tag_ldg.remove_transitive_links(10);

        tag_ldg.report_connectivity();
        std::cout<<"Experimental: linearize bubbles"<<std::endl;
        std::vector<sgNodeID_t> to_remove;
        auto all_paired_linkage = lu.make_paired_linkage(2);
        uint64_t b12=0,b21=0, amb=0, none=0,ab12=0,ab21=0;
        for (auto bubble:tag_ldg.find_bubbles(0,20000)) {
            //if (full_tag_ldg.are_connected(tag_ldg.get_bw_links(bubble.first)[0].dest,tag_ldg.get_fw_links(bubble.first)[0].dest))
            //    to_remove.emplace_back(llabs(bubble.first));
            bool c12=all_paired_linkage.are_connected(-bubble.first,bubble.second);
            bool c21=all_paired_linkage.are_connected(-bubble.second,bubble.first);
            bool i1=all_paired_linkage.are_connected(-bubble.second,-bubble.first);
            bool i2=all_paired_linkage.are_connected(bubble.second,bubble.first);
            if (i1 or i2) {
                ++amb;
                if (c12 and !c21 ) ++ab12;
                else if (!c12 and c21 ) ++ab21;
                continue;
            }
            else if (c12 and !c21 ) {
                tag_ldg.add_link(-bubble.first,bubble.second,0);
                ++b12;
            }
            else if (!c12 and c21 ) {
                tag_ldg.add_link(-bubble.first,bubble.second,0);
                ++b21;
            }
            else {
                ++none;
            }

        }
        std::cout<<"Bubbles with possible reconnection: 1-2: "<<b12<<"  2-1:"<<b21<<" ambiguous order: "<<amb<<"( 1-2: "<<ab12<<" 2-1: "<<ab21<<" ) no connection:"<<none<<std::endl;
        tag_ldg.remove_transitive_links(10);
        tag_ldg.report_connectivity();
        std::cout<<"TODO: clip tips, can this be done??"<<std::endl;
        std::cout<<"TODO: evaluate 'complex nodes' that have a linear transitivity"<<std::endl;
        std::cout<<"TODO: maybe even try local assemblies to simplify complex connections?"<<std::endl;
        ws.sg.write_to_gfa(output_prefix+"_linkage.gfa", {}, {}, selnodes, tag_ldg.links);
        tag_ldg.dump_to_text(dev_create_linkage);
        exit(0);
    }
    if (!dev_skate_linkage.empty()) {
        LinkageUntangler lu(ws);
        LinkageDiGraph tag_ldg(ws.sg);
        tag_ldg.load_from_text(dev_skate_linkage);
        tag_ldg.report_connectivity();
        lu.select_nodes_by_size_and_ci(min_backbone_node_size,min_backbone_ci,max_backbone_ci);
        lu.expand_linear_regions_skating(tag_ldg,dev_max_lines);
        exit(0);
    }
    if (!dev_local_assembly_linkage.empty()) {
        sglib::OutputLog()<<"STARTING DEVEL LOCAL ASSEMBLY RUN"<<std::endl;
        LinkageUntangler lu(ws);
        LinkageDiGraph tag_ldg(ws.sg);
        sglib::OutputLog()<<"Loading linkage from text"<<std::endl;
        tag_ldg.load_from_text(dev_local_assembly_linkage);
        sglib::OutputLog()<<"Analysing connectivity"<<std::endl;
        tag_ldg.report_connectivity();
        sglib::OutputLog()<<"Creating local assembly problems..."<<std::endl;
//        lu.linear_regions_tag_local_assembly(tag_ldg, dev_local_k, dev_local_min_cvg, dev_max_lines,dev_min_nodes,dev_min_total_size,true);
//        ws.sg.write_to_gfa(output_prefix+"_local_patched.gfa");
//        ws.dump_to_disk(output_prefix+"_local_patched.bsgws");
        auto lines=tag_ldg.get_all_lines(dev_min_nodes);
        if (dev_max_lines) lines.resize(dev_max_lines);
        uint64_t li=0;
        uint64_t i=0;
        LinkageUntangler lu2(ws);
        for (auto l:lines) {
            ++i;
            std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
            LocalHaplotypeAssembler lha(ws);
            lha.init_from_backbone(l);
            for (auto ln:l) lu2.selected_nodes[llabs(ln)]=true;
            //lha.assemble(63,5,false);
            //lha.write_problem("local_hap_problem_"+std::to_string(++li));
            //lha.write_full("local_hap_full_"+std::to_string(li));
            lha.assemble(63,5,false);
            lha.write_anchors("local_"+std::to_string(i)+"_anchors.fasta");
            lha.write_gfa("local_"+std::to_string(i)+".gfa");
            lha.construct_patches();
            lha.write_patches("local_"+std::to_string(i)+"_patches.fasta");
            for (auto &p:lha.patches) {

            }
            std::chrono::steady_clock::time_point end= std::chrono::steady_clock::now();
            std::cout << "Local assembly done in " << std::chrono::duration_cast<std::chrono::seconds>(end - begin).count() <<" seconds, produced "<<lha.patches.size()<<" patches"<<std::endl;

        }
        sglib::OutputLog()<<"---NODES CONNECTED ON GLOBAL PROBLEM: "<<std::endl;
        for (auto n=1;n<ws.sg.nodes.size();++n) {
            if (tag_ldg.get_fw_links(n).size()>0 or tag_ldg.get_bw_links(n).size()>0) lu.selected_nodes[n]=true;
        }
        lu.report_node_selection();
        sglib::OutputLog()<<"Bubbles in the linkage digraph: "<<tag_ldg.find_bubbles(0,10000000).size()<<std::endl;
        sglib::OutputLog()<<"---NODES USED ON LOCAL PROBLEMS ("<<lines.size()<<" lines): "<<std::endl;
        lu2.report_node_selection();
        exit(0);
    }
    if (!dev_linkage_stats.empty()) {
        sglib::OutputLog()<<"STARTING DEVEL LOCAL ASSEMBLY RUN"<<std::endl;
        LinkageUntangler lu(ws);
        LinkageDiGraph tag_ldg(ws.sg);
        sglib::OutputLog()<<"Loading linkage from text"<<std::endl;
        tag_ldg.load_from_text(dev_linkage_stats);
        sglib::OutputLog()<<"Analysing connectivity"<<std::endl;
        tag_ldg.report_connectivity();
        sglib::OutputLog()<<"Creating local assembly problems..."<<std::endl;
//        lu.linear_regions_tag_local_assembly(tag_ldg, dev_local_k, dev_local_min_cvg, dev_max_lines,dev_min_nodes,dev_min_total_size,true);
//        ws.sg.write_to_gfa(output_prefix+"_local_patched.gfa");
//        ws.dump_to_disk(output_prefix+"_local_patched.bsgws");
        auto lines=tag_ldg.get_all_lines(dev_min_nodes);
        if (dev_max_lines) lines.resize(dev_max_lines);
        uint64_t li=0;
        LinkageUntangler lu2(ws);
        for (auto l:lines) {
            //LocalHaplotypeAssembler lha(ws);
            //lha.init_from_backbone(l);
            for (auto ln:l) lu2.selected_nodes[llabs(ln)]=true;
            //lha.assemble(63,5,false);
            //lha.write_problem("local_hap_problem_"+std::to_string(++li));
            //lha.write_full("local_hap_full_"+std::to_string(li));
        }
        sglib::OutputLog()<<"---NODES CONNECTED ON GLOBAL PROBLEM: "<<std::endl;
        for (auto n=1;n<ws.sg.nodes.size();++n) {
            if (tag_ldg.get_fw_links(n).size()>0 or tag_ldg.get_bw_links(n).size()>0) lu.selected_nodes[n]=true;
        }
        lu.report_node_selection();
        sglib::OutputLog()<<"Bubbles in the linkage digraph: "<<tag_ldg.find_bubbles(0,10000000).size()<<std::endl;
        sglib::OutputLog()<<"---NODES USED ON LOCAL PROBLEMS ("<<lines.size()<<" lines): "<<std::endl;
        lu2.report_node_selection();
        exit(0);
    }
    if (dev_local_patching) {
        LinkageUntangler lu(ws);
        ws.kci.reindex_graph();
        if (select_hspnps) lu.select_nodes_by_HSPNPs(min_backbone_node_size,min_backbone_ci,max_backbone_ci);
        else lu.select_nodes_by_size_and_ci(min_backbone_node_size,min_backbone_ci,max_backbone_ci);
        std::unordered_set<sgNodeID_t> selnodes;
        for (sgNodeID_t n=1;n<ws.sg.nodes.size();++n) if (lu.selected_nodes[n]) selnodes.insert(n);
        lu.report_node_selection();
        auto pre_tag_ldg = lu.make_tag_linkage(min_shared_tags,dev_linkage_paths);
        pre_tag_ldg.remove_transitive_links(10);
        pre_tag_ldg.report_connectivity();
        sglib::OutputLog()<<"Eliminating N-N nodes..."<<std::endl;
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
        ws.sg.write_to_gfa(output_prefix+"_linkage.gfa", {}, {}, selnodes, tag_ldg.links);
        lu.linear_regions_tag_local_assembly(tag_ldg, dev_local_k, dev_local_min_cvg, dev_max_lines,dev_min_nodes,dev_min_total_size,true);
    }

    ws.kci.reindex_graph();
    if (dump_gfa) ws.sg.write_to_gfa(output_prefix+".gfa");
    if (remap_reads) ws.remap_all();
    ws.dump_to_disk(output_prefix+".bsgws");
    return 0;
}

