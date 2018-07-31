#include <iostream>
#include <fstream>
#include <sglib/WorkSpace.hpp>
#include <sglib/processors/Untangler.hpp>
#include <sglib/processors/FlowFollower.hpp>
#include <sglib/processors/LinkageUntangler.hpp>
#include <sglib/processors/LocalHaplotypeAssembler.hpp>
#include <sglib/processors/GraphEditor.hpp>
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
    std::string create_linkage,dev_skate_linkage,make_patches, load_patches, patch_backbones, patch_workspace,load_linkage;
    bool dev_local_patching=false;
    bool remap_reads=true;
    bool dump_gfa=false;
    bool dev_linkage_paths=false;
    bool dev_test_make_patches=false;
    std::string dev_linkage_stats;
    int dev_dump_local_problems_from=-1;
    int dev_dump_local_problems_to=-1;
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

        options.add_options("Untangling")
                ("p,paired_reads","paired read scaffolding (experimental)",cxxopts::value<bool>(paired_scaff))
                ("l,unroll_loops", "unroll simple loops",cxxopts::value<bool>(unroll_loops))
                ("e,pop_errors", "pop unsupported short-bubbles (as errors)",cxxopts::value<bool>(pop_errors))
                ("r,repeat_expansion","run tag-based repeat expansion", cxxopts::value<bool>(repeat_expansion))
                ("b,bubbly_paths","run bubbly paths phasing", cxxopts::value<bool>(bubbly_paths))
                ("min_pairs","minimum number of pairs to connect two nodes",cxxopts::value<int>(min_pairs))
                ("tag_imbalance_ends","percentage of node to use as tag-imbalanced end",cxxopts::value<float>(tag_imbalance_ends))
                ("min_shared_tags","minimum shared tags to evaluate tag-imbalanced connection",cxxopts::value<int>(min_shared_tags));

        options.add_options("Linkage untangling")
                ("create_linkage","Creates and simplifies linkage and dumps to file",cxxopts::value<std::string>(create_linkage))
                ("load_linkage","Linkage file to load",cxxopts::value<std::string>(load_linkage))
                ("make_patches","Creates patches from local assemblies",cxxopts::value<std::string>(make_patches))
                ("load_patches","Load patches form file",cxxopts::value<std::string>(load_patches))
                ("patch_backbones","Patches backbones and outputs stitched sequences",cxxopts::value<std::string>(patch_backbones))
                ("patch_workspace","Patches  the workspace graph and outputs a gfa",cxxopts::value<std::string>(patch_workspace))
                ("full_local_patching","run a whole round of linked lines and local patching, with a final remap", cxxopts::value<bool>(dev_local_patching));
        options.add_options("Development")

                ("dev_linkage_paths", "tag linkage uses read pathing rather than simple mapping",cxxopts::value<bool>(dev_linkage_paths))
                //("dev_skate_linkage","Loads linkage from file and skates",cxxopts::value<std::string>(dev_skate_linkage))
                ("dev_linkage_stats","Loads linkage from file and computes local assemblies stats",cxxopts::value<std::string>(dev_linkage_stats))
                ("dev_max_lines","Limits lines to be skated on dev",cxxopts::value<int>(dev_max_lines))
                ("dev_min_nodes","Limits lines to be locally assembled on dev to at least min_nodes",cxxopts::value<uint64_t>(dev_min_nodes))
                ("dev_min_total_size","Limits lines to be locally assembled on dev to at least min_total_size",cxxopts::value<uint64_t>(dev_min_total_size))
                ("dev_local_k","k value for local assembly",cxxopts::value<uint8_t>(dev_local_k))
                ("dev_local_min_cvg","minimum coverga for local assemblies",cxxopts::value<int>(dev_local_min_cvg))
                ("dev_dump_local_problems_from","dumps local problems from",cxxopts::value<int>(dev_dump_local_problems_from))
                ("dev_dump_local_problems_to","dumps local problems from",cxxopts::value<int>(dev_dump_local_problems_to));




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

    //======= DEVELOPMENT FUNCTIONS ======
    if (!create_linkage.empty()) {
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
        tag_ldg.dump_to_text(create_linkage);
        exit(0);
    }

    LinkageDiGraph tag_ldg(ws.sg);
    if (!load_linkage.empty()){
        tag_ldg.load_from_text(load_linkage);
    }

    if (dev_dump_local_problems_from>-1 and dev_dump_local_problems_from<dev_dump_local_problems_to){
        sglib::OutputLog()<<"Analysing connectivity"<<std::endl;
        tag_ldg.report_connectivity();
        sglib::OutputLog()<<"Creating local assembly problems..."<<std::endl;
        auto lines=tag_ldg.get_all_lines(dev_min_nodes);
        if (dev_max_lines) {
            sglib::OutputLog()<<"dev_max_lines set, solving only "<<dev_max_lines<<" / "<<lines.size()<<std::endl;
            lines.resize(dev_max_lines);
        }
        uint64_t li=0;
        uint64_t i=0;
        sglib::OutputLog()<<"Dumping from "<< dev_dump_local_problems_from <<" to "<< dev_dump_local_problems_to << " of " <<lines.size()<<" local assembly problems..."<<std::endl;
        for (auto li=dev_dump_local_problems_from;li<lines.size() and li<=dev_dump_local_problems_to;++li) {
            auto &l = lines[li];
            LocalHaplotypeAssembler lha(ws);
            lha.init_from_backbone(l);
            lha.write_full(output_prefix + "_local_" + std::to_string(li) + ".bsglhapf");
        }
        exit(0);
    }

    std::vector<std::pair<std::pair<sgNodeID_t ,sgNodeID_t>,std::string>> patches;
    if (!make_patches.empty()) {
        sglib::OutputLog()<<"Analysing connectivity"<<std::endl;
        tag_ldg.report_connectivity();
        sglib::OutputLog()<<"Creating local assembly problems..."<<std::endl;
        auto lines=tag_ldg.get_all_lines(dev_min_nodes);
        if (dev_max_lines) {
            sglib::OutputLog()<<"dev_max_lines set, solving only "<<dev_max_lines<<" / "<<lines.size()<<std::endl;
            lines.resize(dev_max_lines);
        }
        uint64_t li=0;
        uint64_t i=0;
        sglib::OutputLog()<<"Solving "<<lines.size()<<" local assembly problems..."<<std::endl;
#pragma omp parallel for schedule(dynamic,10)
        for (auto li=0;li<lines.size();++li) {
            auto &l=lines[li];
            std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
            LocalHaplotypeAssembler lha(ws);
            lha.init_from_backbone(l);
            lha.assemble(63,5,false,false);
            lha.construct_patches();
            //lha.write_patches("local_"+std::to_string(i)+"_patches.fasta");
            std::chrono::steady_clock::time_point end= std::chrono::steady_clock::now();
#pragma omp critical(patch_collection)
            {
                patches.insert(patches.end(),lha.patches.begin(),lha.patches.end());
                sglib::OutputLog() << "Local assembly #"<<li<<" from " << l.size() << " anchors done in "
                                   << std::chrono::duration_cast<std::chrono::seconds>(end - begin).count()
                                   << " seconds, produced " << lha.patches.size() << " patches" << std::endl;
            }
        }

        sglib::OutputLog()<<patches.size()<<" patches found"<<std::endl;
        std::ofstream patchf(make_patches);
        for (auto &p:patches){
            patchf << ">patch_" << -p.first.first << "_" << p.first.second << std::endl;
            patchf << p.second << std::endl;
        }
        exit(0);
    }

    if (!load_patches.empty()){
        std::ifstream pf(load_patches);
        auto i=0;
        while (!pf.eof()){
            std::string l;
            pf>>l;
            if (l.empty()) break;
            l=l.substr(7,l.size());
            auto l2=l.substr(l.find('_')+1,l.size());
            sgNodeID_t n1=atol(l.c_str());
            sgNodeID_t n2=atol(l2.c_str());
            std::string seq;
            pf>>seq;
            patches.emplace_back(std::make_pair(n1,n2),seq);
        }
    }

    if (!dev_linkage_stats.empty()) {
        LinkageUntangler lu(ws);
        sglib::OutputLog()<<"Analysing connectivity"<<std::endl;
        tag_ldg.report_connectivity();
        auto lines=tag_ldg.get_all_lines(dev_min_nodes);
        if (dev_max_lines) lines.resize(dev_max_lines);
        uint64_t li=0;
        LinkageUntangler lu2(ws);
        for (auto l:lines) {
            for (auto ln:l) lu2.selected_nodes[llabs(ln)]=true;
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

    if (!patch_backbones.empty()) {
        std::ofstream patched_backbones_file(patch_backbones);
        sglib::OutputLog() << "Analysing connectivity" << std::endl;
        tag_ldg.report_connectivity();
        sglib::OutputLog() << "Creating local assembly problems..." << std::endl;
        auto lines = tag_ldg.get_all_lines(dev_min_nodes);
        if (dev_max_lines) {
            sglib::OutputLog() << "dev_max_lines set, solving only " << dev_max_lines << " / " << lines.size()
                               << std::endl;
            lines.resize(dev_max_lines);
        }
        uint64_t li = 0;
        uint64_t i = 0;
        std::map<std::pair<sgNodeID_t, sgNodeID_t>, std::string> patchmap;
        for (auto p:patches) patchmap[p.first]=p.second;
        sglib::OutputLog() << "Patching in " << lines.size() << " local assembly problems from " << patches.size() << " patches..." << std::endl;
        int endsize=200;
        for (auto li = 0; li < lines.size(); ++li) {
            auto &l = lines[li];
            std::string name=">"+std::to_string(l[0]);
            auto n=ws.sg.nodes[llabs(l[0])];
            if (l[0]<0)n.make_rc();
            std::string seq=n.sequence;
            for (auto ni=0;ni<l.size()-1;++ni){
                //for each transition, either put the patch, or start a new sequence
                if (patchmap.count(std::make_pair(-l[ni],l[ni+1]))==0){
                    //patch not found
                    patched_backbones_file<<name<<std::endl<<seq<<std::endl;
                    name=">"+std::to_string(l[ni+1]);
                    n=ws.sg.nodes[llabs(l[ni+1])];
                    if (l[ni+1]<0)n.make_rc();
                    seq=n.sequence;
                } else {

                    //patch found
                    //create/add patch substring
                    std::string p=patchmap[std::make_pair(-l[ni],l[ni+1])];
                    seq+=p.substr(endsize+1,p.size()-endsize*3-1);
                    //add next node
                    name+="_"+std::to_string(l[ni+1]);
                    n=ws.sg.nodes[llabs(l[ni+1])];
                    if (l[ni+1]<0)n.make_rc();
                    seq+=n.sequence;
                }
            }
            patched_backbones_file<<name<<std::endl<<seq<<std::endl;
        }
        exit(0);
    }

    if (!patch_workspace.empty()) {
        sglib::OutputLog() << "Patching workspace!!!" << std::endl;
        GraphEditor ge(ws);
        uint64_t patch_results[6]={0,0,0,0,0,0};
        for (auto p:patches) {
            auto r=ge.patch_between(-p.first.first,p.first.second,p.second);
            ++patch_results[r];
        }
        sglib::OutputLog() << "Patches with no anchor ends:         "<< patch_results[1] <<std::endl;
        sglib::OutputLog() << "Patches with incorrect anchor order: "<< patch_results[3] <<std::endl;
        sglib::OutputLog() << "Patches with no SG paths:            "<< patch_results[2] <<std::endl;
        sglib::OutputLog() << "Patches with no matching SG paths:   "<< patch_results[4] <<std::endl;
        sglib::OutputLog() << "Patches with failed expansion:       "<< patch_results[5] <<std::endl;
        sglib::OutputLog() << "Patches correctly applied:           "<< patch_results[0] <<std::endl;
        auto juc=ws.sg.join_all_unitigs();
        sglib::OutputLog() << juc << " unitigs joined after patching"<<std::endl;
        ws.sg.write_to_gfa(patch_workspace);
    }


    ws.kci.reindex_graph();
    if (dump_gfa) ws.sg.write_to_gfa(output_prefix+".gfa");
    if (remap_reads) ws.remap_all();
    ws.dump_to_disk(output_prefix+".bsgws");
    return 0;
}

