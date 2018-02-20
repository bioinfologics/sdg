#include <iostream>
#include <fstream>
#include <sglib/PairedReadMapper.h>
#include <sglib/Scaffolder.hpp>
#include <sglib/KmerCompressionIndex.hpp>
#include <sglib/GraphPartitioner.hpp>
#include <sglib/mappers/LinkedReadMapper.hpp>
#include "sglib/SequenceGraph.h"
#include "sglib/logger/OutputLog.h"
#include "cxxopts.hpp"


int main(int argc, char * argv[]) {
    std::cout << "Welcome to gfa-pairscaff"<<std::endl<<std::endl;
    std::cout << "Git origin: " << GIT_ORIGIN_URL << " -> "  << GIT_BRANCH << std::endl;
    std::cout << "Git commit: " << GIT_COMMIT_HASH << std::endl<<std::endl;
    std::cout << "Executed command:"<<std::endl;
    for (auto i=0;i<argc;i++) std::cout<<argv[i]<<" ";
    std::cout<<std::endl<<std::endl;

    std::string gfa_filename,output_prefix, load_cidx, dump_cidx;
    std::vector<std::string> reads1,reads2,reads_type,cidxreads1,cidxreads2, dump_mapped, load_mapped;
    bool stats_only=0,reindex_kci=0;
    uint64_t max_mem_gb=4;
    sglib::OutputLogLevel=sglib::LogLevels::DEBUG;
    bool repeats_first=1;
    try
    {
        cxxopts::Options options("gfa-bfphaser", "GFA Bubbles First Phaser");

        options.add_options()
                ("help", "Print help")
                ("g,gfa", "input gfa file", cxxopts::value<std::string>(gfa_filename))
                ("o,output", "output file prefix", cxxopts::value<std::string>(output_prefix));
        options.add_options("Compression index options")
                ("cidxread1", "compression index input reads, left", cxxopts::value<std::vector<std::string>>(cidxreads1))
                ("cidxread2", "compression index input reads, right", cxxopts::value<std::vector<std::string>>(cidxreads2))
                ("load_cidx", "load compression index filename", cxxopts::value<std::string>(load_cidx))
                ("dump_cidx", "dump compression index filename", cxxopts::value<std::string>(dump_cidx))
                ("reindex_ci", "re-index compression index after loading", cxxopts::value<bool>(reindex_kci))
                ;
        options.add_options("Paired reads options")
                ("1,read1", "input reads, left", cxxopts::value<std::vector<std::string>>(reads1))
                ("2,read2", "input reads, right", cxxopts::value<std::vector<std::string>>(reads2))
                ("read_type", "One of: pe,10x", cxxopts::value<std::vector<std::string>>(reads_type))
                ("d,dump_to", "dump mapped reads to file", cxxopts::value<std::vector<std::string>>(dump_mapped))
                ("l,load_from", "load mapped reads from file", cxxopts::value<std::vector<std::string>>(load_mapped))
                ("max_mem", "maximum_memory when mapping (GB, default: 4)", cxxopts::value<uint64_t>(max_mem_gb));

        options.add_options("Phasing and scaffolding options")
                ("repeats_first", "solve_repeats before bubble phasing", cxxopts::value<bool>(repeats_first))
                ;



        auto result(options.parse(argc, argv));

        if (result.count("help"))
        {
            std::cout << options.help({"","Paired reads options","Compression index options","Phasing and scaffolding options"}) << std::endl;
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

    std::cout<<std::endl<<"=== Loading GFA ==="<<std::endl;
    if (gfa_filename.size()<=4 or gfa_filename.substr(gfa_filename.size()-4,4)!=".gfa") {

        throw std::invalid_argument("filename of the gfa input does not end in gfa, it ends in '" +
                                    gfa_filename.substr(gfa_filename.size() - 4, 4) + "'");
    }
    auto fasta_filename=gfa_filename.substr(0,gfa_filename.size()-4)+".fasta";
    SequenceGraph sg;
    sg.load_from_gfa(gfa_filename);

    std::cout<<std::endl<<"=== Mapping reads ==="<<std::endl;
    //read mapping/loading
    std::vector<PairedReadMapper> prmappers;
    std::vector<LinkedReadMapper> lrmappers;
    std::vector<LinkedReadsDatastore> datastores;


    for(int lib=0;lib<reads1.size();lib++) {
        if (load_mapped.size()<=lib) {
            if (reads_type[lib] == "10x" or reads_type[lib] == "10xseq") {
                datastores.emplace_back(reads1[lib],reads2[lib],
                        (reads_type[lib] == "10xseq" ? LinkedReadsFormat::seq : LinkedReadsFormat::UCDavis));
                lrmappers.emplace_back(sg,datastores.back());
                lrmappers.back().memlimit=max_mem_gb * 1024L * 1024L * 1024L;
                lrmappers.back().update_graph_index();
                lrmappers.back().map_reads();
                //mappers.back().map_reads(reads1[lib], reads2[lib], prm10x, max_mem_gb * 1024L * 1024L * 1024L);
            } else {
                std::cout<<"Sorry, but only 10x read mapping is supported at this time"<<std::endl;
            }
            lrmappers.back().print_stats();
        } else {
            std::cout<<"Sorry, loading datastores and mappers is NOT supported at this time"<<std::endl;
            /*std::cout<<"Library #"<<lib<<": load from disk at "<<load_mapped[lib]<<std::endl;
            mappers.emplace_back(sg);
            mappers.back().load_from_disk(load_mapped[lib]);
            mappers[lib].read1filename=reads1[lib];
            mappers[lib].read2filename=reads2[lib];
            mappers[lib].readType=(reads_type[lib] == "10x" ? prm10x : prmPE);
            mappers[lib].memlimit=max_mem_gb * 1024L * 1024L * 1024L;
            mappers.back().print_stats();*/
        }
        if (dump_mapped.size() > lib) {
            std::cout<<"Sorry, dumping datastores and mappers is NOT supported at this time"<<std::endl;
            /*std::cout<<"dumping map to "<<dump_mapped[lib]<<std::endl;
            mappers.back().save_to_disk(dump_mapped[lib]);*/
        }
    }

    std::cout<<std::endl<<"=== Loading reads compression index ==="<<std::endl;
    //compression index
    KmerCompressionIndex kci(sg,max_mem_gb*1024L*1024L*1024L);
    if (load_cidx!=""){
        kci.load_from_disk(load_cidx);
        if (reindex_kci) kci.reindex_graph();
    } else {
        kci.index_graph();
        for(int lib=0;lib<cidxreads1.size();lib++) {
            kci.start_new_count();
            kci.add_counts_from_file(cidxreads1[lib]);
            kci.add_counts_from_file(cidxreads2[lib]);
        }
    }
    if (dump_cidx!=""){
        kci.save_to_disk(dump_cidx);
    }

    if (kci.read_counts.size()>0) {
        kci.compute_compression_stats();
        kci.dump_histogram("kci_histogram.csv");
    }


    std::cout<<std::endl<<"=== Scaffolding ==="<<std::endl;


    Scaffolder scaff(sg,prmappers,lrmappers,kci);

    if (repeats_first) {
        std::cout << std::endl << "Step 1 - Solving trivial repeats" << std::endl;
        bool mod = true;
        int srpass = 0;
        while (mod) {
            uint64_t aa_count=0,ab_count=0,unsolved_count=0;
            mod = false;
            ++srpass;
            std::cout << " Finding trivial repeats to analyse with tags" << std::endl;
            std::vector<bool> used(sg.nodes.size());
            std::vector<SequenceGraphPath> paths_solved;
            for (auto n = 1; n < sg.nodes.size(); ++n) {
                if (used[n]) continue;
                auto fwl = sg.get_fw_links(n);
                auto bwl = sg.get_bw_links(n);
                if (fwl.size() != 2 or bwl.size() != 2) continue;
                auto f0 = fwl[0].dest;
                auto f1 = fwl[1].dest;
                auto b0 = bwl[0].dest;
                auto b1 = bwl[1].dest;
                if (used[(b0 > 0 ? b0 : -b0)] or used[(b1 > 0 ? b1 : -b1)] or used[(f0 > 0 ? f0 : -f0)] or
                    used[(f1 > 0 ? f1 : -f1)])
                    continue;
                sgNodeID_t all[4] = {f0, f1, b0, b1};
                bool ok = true;
                for (auto x:all) if (x == n or x == -n or sg.nodes[(x > 0 ? x : -x)].sequence.size() < 399) ok = false;
                for (auto j = 0; j < 3; ++j)
                    for (auto i = j + 1; i < 4; ++i)
                        if (all[i] == all[j] or all[i] == -all[j])ok = false; //looping node
                if (!ok) continue;

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
                    paths_solved.push_back(SequenceGraphPath(sg, {-b0, n, f0}));
                    paths_solved.push_back(SequenceGraphPath(sg, {-b1, n, f1}));
                    ++aa_count;

                } else if (ba.size() > 3 and ab.size() > 3 and
                           std::min(ba.size(), ab.size()) > 10 * std::max(aa.size(), bb.size())) {
                    //std::cout << " Solved as AB BA !!!" << std::endl;
                    used[(b0 > 0 ? b0 : -b0)] = true;
                    used[(b1 > 0 ? b1 : -b1)] = true;
                    used[n] = true;
                    used[(f0 > 0 ? f0 : -f0)] = true;
                    used[(f1 > 0 ? f1 : -f1)] = true;
                    paths_solved.push_back(SequenceGraphPath(sg, {-b0, n, f1}));
                    paths_solved.push_back(SequenceGraphPath(sg, {-b1, n, f0}));
                    ++ab_count;
                }
                else ++unsolved_count;
            }

            std::cout << paths_solved.size() << " paths to join" << std::endl;
            if (paths_solved.size() > 0) mod = true;
            std::unordered_set<uint64_t> reads_to_remap;
            for (auto &p:paths_solved) {
                sg.join_path(p);
                for (auto n:p.nodes) {
                    for (auto rm:scaff.lrmappers[0].reads_in_node[(n > 0 ? n : -n)]) {
                        reads_to_remap.insert((rm.read_id % 2 ? rm.read_id : rm.read_id - 1));
                        scaff.lrmappers[0].read_to_node[rm.read_id] = 0;
                    }
                }

            }
            std::cout<<"Path analysis summary AA:"<<aa_count<<" AB:"<<ab_count<<" Unsolved:"<<unsolved_count<<std::endl;
            if (mod) {
                scaff.kci.reindex_graph();
                sg.write_to_gfa(output_prefix + "_solved_repeats_" + std::to_string(srpass) + ".gfa");
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






    sg.write_to_gfa(output_prefix+"_scaffolded.gfa");
    return 0;
}

