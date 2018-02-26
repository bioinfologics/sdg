#include <iostream>
#include <fstream>
#include <sglib/WorkSpace.hpp>
#include <sglib/GraphPartitioner.hpp>
#include <sglib/Scaffolder.hpp>
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
    bool repeats_first=1,pop_bubles=0,skip_remap=0;
    std::string verbose_log="";
    try
    {
        cxxopts::Options options("bsg-untangler", "graph-based repeat resolution and haplotype separation");

        options.add_options()
                ("help", "Print help")
                ("w,workspace", "input workspace", cxxopts::value<std::string>(workspace_file))
                ("o,output", "output file prefix", cxxopts::value<std::string>(output_prefix));
                ("max_mem", "maximum_memory when mapping (GB, default: 4)", cxxopts::value<uint64_t>(max_mem_gb));

        options.add_options("Heuristics")
                ("skip_remap", "skip remapping the reads after loading workspace", cxxopts::value<bool>(skip_remap))
                ("pop_short", "pop unsupported short bubbles first", cxxopts::value<bool>(pop_bubles))
                ("repeats_first", "solve_repeats before bubble phasing", cxxopts::value<bool>(repeats_first))
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


    WorkSpace ws;
    std::cout<<std::endl<<"=== Loading Workspace ==="<<std::endl;
    ws.load_from_disk(workspace_file);

    std::cout<<std::endl<<"=== Mapping reads ==="<<std::endl;
    ws.linked_read_mappers[0].memlimit=max_mem_gb*1024*1024*1024;
    if (skip_remap==false) {
        ws.linked_read_mappers[0].update_graph_index();
        ws.linked_read_mappers[0].map_reads();
        ws.dump_to_disk(output_prefix+"_mapped.bsgws");
    }




    std::cout<<std::endl<<"=== Scaffolding ==="<<std::endl;

    std::vector<PairedReadMapper> prmappers;
    Scaffolder scaff(ws.sg,prmappers,ws.linked_read_mappers,ws.kci);
    std::ofstream verbose_log_file;
    if (verbose_log!=""){
        verbose_log_file.open(verbose_log);
    }
    if (pop_bubles) {
        scaff.pop_unsupported_shortbubbles();
        //for (auto & rm: scaff.rmappers) {
        //    rm.remove_obsolete_mappings();
        //    rm.update_graph_index();
        //    rm.map_reads();
        //}
        for (auto & lrm: scaff.lrmappers) {
            lrm.remove_obsolete_mappings();
            lrm.update_graph_index();
            lrm.map_reads();
        }
    }
    std::unordered_set<sgNodeID_t> unsolved_repeats;
    if (repeats_first) {
        std::cout << std::endl << "Step 1 - Solving trivial repeats" << std::endl;
        bool mod = true;
        int srpass = 0;
        while (mod) {
            unsolved_repeats.clear();
            uint64_t aa_count=0,ab_count=0,unsolved_count=0,non_evaluated=0;
            mod = false;
            ++srpass;
            std::cout << " Finding trivial repeats to analyse with tags" << std::endl;
            std::vector<bool> used(ws.sg.nodes.size());
            std::vector<SequenceGraphPath> paths_solved;
            if (verbose_log!="") verbose_log_file<<"==== Round of repetition analysis started ===="<<std::endl;
            for (auto n = 1; n < ws.sg.nodes.size(); ++n) {
                if (used[n]) {
                    if (verbose_log!="") verbose_log_file<<n<<" is used"<<std::endl;
                    continue;
                }
                auto fwl = ws.sg.get_fw_links(n);
                auto bwl = ws.sg.get_bw_links(n);
                if (fwl.size() != 2 or bwl.size() != 2) {
                    if (verbose_log!="") verbose_log_file<<n<<" has "<<bwl.size()<<" ins and "<<fwl.size()<<" outs"<<std::endl;
                    continue;
                }
                auto f0 = fwl[0].dest;
                auto f1 = fwl[1].dest;
                auto b0 = bwl[0].dest;
                auto b1 = bwl[1].dest;
                if (used[(b0 > 0 ? b0 : -b0)] or used[(b1 > 0 ? b1 : -b1)] or used[(f0 > 0 ? f0 : -f0)] or
                    used[(f1 > 0 ? f1 : -f1)]) {
                    non_evaluated++;
                    if (verbose_log!="") verbose_log_file<<n<<" neighbors are used"<<std::endl;
                    continue;
                }
                sgNodeID_t all[4] = {f0, f1, b0, b1};
                bool ok = true;
                for (auto x:all) if (x == n or x == -n or ws.sg.nodes[(x > 0 ? x : -x)].sequence.size() < 199) ok = false;
                for (auto j = 0; j < 3; ++j)
                    for (auto i = j + 1; i < 4; ++i)
                        if (all[i] == all[j] or all[i] == -all[j])ok = false; //looping node
                if (!ok) {
                    if (verbose_log!="") verbose_log_file<<n<<" and neighbors are short"<<std::endl;
                    continue;
                }

//                std::cout<<"Repeat: [ "<<-b0<<" | "<<-b1<<" ] <-> "<<n<<" <-> [ "<<f0<<" | "<<f1<<" ]"<<std::endl;
                std::set<prm10xTag_t> b0tags, b1tags, f0tags, f1tags;

//                std::cout<<"Reads in nodes -> f0:"<<scaff.lrmappers[0].reads_in_node[(f0 > 0 ? f0 : -f0)].size()
//                        <<"   f1:"<<scaff.lrmappers[0].reads_in_node[(f1 > 0 ? f1 : -f1)].size()
//                        <<"   b0:"<<scaff.lrmappers[0].reads_in_node[(b0 > 0 ? b0 : -b0)].size()
//                        <<"   b1:"<<scaff.lrmappers[0].reads_in_node[(b1 > 0 ? b1 : -b1)].size()<<std::endl;

                for (auto rm:scaff.lrmappers[0].reads_in_node[(f0 > 0 ? f0 : -f0)])
                    f0tags.insert(scaff.lrmappers[0].datastore.get_read_tag(rm.read_id));
                for (auto rm:scaff.lrmappers[0].reads_in_node[(f1 > 0 ? f1 : -f1)])
                    f1tags.insert(scaff.lrmappers[0].datastore.get_read_tag(rm.read_id));
                for (auto rm:scaff.lrmappers[0].reads_in_node[(b0 > 0 ? b0 : -b0)])
                    b0tags.insert(scaff.lrmappers[0].datastore.get_read_tag(rm.read_id));
                for (auto rm:scaff.lrmappers[0].reads_in_node[(b1 > 0 ? b1 : -b1)])
                    b1tags.insert(scaff.lrmappers[0].datastore.get_read_tag(rm.read_id));

                f0tags.erase(0);
                f1tags.erase(0);
                b0tags.erase(0);
                b1tags.erase(0);

                std::set<prm10xTag_t> shared1,shared2;

                for (auto t:f0tags) if (f1tags.count(t)>0) {shared1.insert(t);};
                for (auto t:shared1){f0tags.erase(t);f1tags.erase(t);};
                for (auto t:b0tags) if (b1tags.count(t)>0) {shared2.insert(t);};
                for (auto t:shared2) {b0tags.erase(t);b1tags.erase(t);};

                std::set<prm10xTag_t> aa, bb, ba, ab;
                std::set_intersection(b0tags.begin(), b0tags.end(), f0tags.begin(), f0tags.end(),
                                      std::inserter(aa, aa.end()));
                std::set_intersection(b0tags.begin(), b0tags.end(), f1tags.begin(), f1tags.end(),
                                      std::inserter(ab, ab.end()));
                std::set_intersection(b1tags.begin(), b1tags.end(), f0tags.begin(), f0tags.end(),
                                      std::inserter(ba, ba.end()));
                std::set_intersection(b1tags.begin(), b1tags.end(), f1tags.begin(), f1tags.end(),
                                      std::inserter(bb, bb.end()));
//                std::cout<<"Tags in   f0: "<<f0tags.size()<<"  f1: "<<f1tags.size()<<"  b0: "<<b0tags.size()<<"  b1: "<<b1tags.size()<<std::endl;
//                std::cout<<"Tags support   aa: "<<aa.size()<<"  bb: "<<bb.size()<<"  ab: "<<ab.size()<<"  ba: "<<ba.size();
                if (aa.size() > 3 and bb.size() > 3 and
                    std::min(aa.size(), bb.size()) > 10 * std::max(ab.size(), ba.size())) {
                    //std::cout << " Solved as AA BB !!!" << std::endl;
                    used[(b0 > 0 ? b0 : -b0)] = true;
                    used[(b1 > 0 ? b1 : -b1)] = true;
                    used[n] = true;
                    used[(f0 > 0 ? f0 : -f0)] = true;
                    used[(f1 > 0 ? f1 : -f1)] = true;
                    paths_solved.push_back(SequenceGraphPath(ws.sg, {-b0, n, f0}));
                    paths_solved.push_back(SequenceGraphPath(ws.sg, {-b1, n, f1}));
                    ++aa_count;

                } else if (ba.size() > 3 and ab.size() > 3 and
                           std::min(ba.size(), ab.size()) > 10 * std::max(aa.size(), bb.size())) {
                    //std::cout << " Solved as AB BA !!!" << std::endl;
                    used[(b0 > 0 ? b0 : -b0)] = true;
                    used[(b1 > 0 ? b1 : -b1)] = true;
                    used[n] = true;
                    used[(f0 > 0 ? f0 : -f0)] = true;
                    used[(f1 > 0 ? f1 : -f1)] = true;
                    paths_solved.push_back(SequenceGraphPath(ws.sg, {-b0, n, f1}));
                    paths_solved.push_back(SequenceGraphPath(ws.sg, {-b1, n, f0}));
                    ++ab_count;
                }
                else {
                    unsolved_repeats.insert(n);
                    ++unsolved_count;
                    if (verbose_log!="") verbose_log_file<<n<<" unsolved aa: "<<aa.size()<<"  bb: "<<bb.size()<<"  ab: "<<ab.size()<<"  ba: "<<ba.size()<<std::endl;
                }
            }

            std::cout << paths_solved.size() << " paths to join" << std::endl;
            if (paths_solved.size() > 0) mod = true;
            std::unordered_set<uint64_t> reads_to_remap;
            for (auto &p:paths_solved) {
                ws.sg.join_path(p);
                for (auto n:p.nodes) {
                    for (auto rm:scaff.lrmappers[0].reads_in_node[(n > 0 ? n : -n)]) {
                        reads_to_remap.insert((rm.read_id % 2 ? rm.read_id : rm.read_id - 1));
                        scaff.lrmappers[0].read_to_node[rm.read_id] = 0;
                    }
                }

            }
            std::cout<<"Path analysis summary AA:"<<aa_count<<" AB:"<<ab_count
                     <<" Unsolved:"<<unsolved_count<<" Non evaluated:"<<non_evaluated<<std::endl;
            if (mod) {
                scaff.kci.reindex_graph();
                ws.sg.write_to_gfa(output_prefix + "_solved_repeats_" + std::to_string(srpass) + ".gfa",unsolved_repeats);
                std::cout << reads_to_remap.size() << " pairs of read to remap" << std::endl;
                scaff.lrmappers[0].update_graph_index();
                scaff.lrmappers[0].map_reads(reads_to_remap);
            }
        }
    }

    std::cout<<std::endl<<"Step 2 - GraphPartitioner"<<std::endl;

    auto hps=scaff.get_all_haplotype_pairs();
    std::cout<<"Starting with "<<hps.size()<<" haplotype pairs"<<std::endl;
    std::ofstream hpsfile(output_prefix+"_hps.csv");
    auto hspnp_id=1;
    for (auto hp:hps) { //TODO: print node numbers and alternative tags -> then cluster/find neighbours
        std::vector<prm10xTag_t> h1tags,h2tags;
        for(auto rm:scaff.lrmappers[0].reads_in_node[(hp.first>0?hp.first:-hp.first)]) h1tags.push_back(scaff.lrmappers[0].datastore.get_read_tag(rm.read_id));
        for(auto rm:scaff.lrmappers[0].reads_in_node[(hp.second>0?hp.second:-hp.second)]) h2tags.push_back(scaff.lrmappers[0].datastore.get_read_tag(rm.read_id));
        std::sort(h1tags.begin(),h1tags.end());
        std::sort(h2tags.begin(),h2tags.end());
        hpsfile<<hspnp_id<<",A,"<<hp.first;
        for (auto t:h1tags) if (t) hpsfile<<","<<t;
        hpsfile<<std::endl;
        hpsfile<<hspnp_id<<",B,"<<hp.second;
        for (auto t:h2tags) if (t) hpsfile<<","<<t;
        hpsfile<<std::endl;
        ++hspnp_id;
    }





    ws.dump_to_disk(output_prefix+"_final.bsgws");
    ws.sg.write_to_gfa(output_prefix+"_scaffolded.gfa");
    return 0;
}

