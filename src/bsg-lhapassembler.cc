#include <iostream>
#include <fstream>
#include <sglib/workspace/WorkSpace.hpp>
#include <sglib/processors/LocalHaplotypeAssembler.hpp>
#include "sglib/logger/OutputLog.hpp"
#include "cxxopts.hpp"


int main(int argc, char * argv[]) {
    std::cout << "Welcome to bsg-lhapassembler"<<std::endl<<std::endl;
    std::cout << "Git origin: " << GIT_ORIGIN_URL << " -> "  << GIT_BRANCH << std::endl;
    std::cout << "Git commit: " << GIT_COMMIT_HASH << std::endl<<std::endl;
    std::cout << "Executed command:"<<std::endl;
    for (auto i=0;i<argc;i++) std::cout<<argv[i]<<" ";
    std::cout<<std::endl<<std::endl;

    sglib::OutputLogLevel=sglib::LogLevels::DEBUG;

    std::string workspace_file,problem_file,output_prefix,benchmark_filename,problem_analysis;

    uint8_t k=63;
    int min_cvg=5;
    bool use_tag_cvg=false;
    try
    {
        cxxopts::Options options("bsg-lhapassembler", "Local Haplotype(specific) Assembler");

        options.add_options()
                ("help", "Print help")
                ("w,workspace", "input workspace", cxxopts::value<std::string>(workspace_file))
                ("p,problem", "problem file", cxxopts::value<std::string>(problem_file))
                ("o,output", "output file prefix", cxxopts::value<std::string>(output_prefix));

        options.add_options("Assembly")
                ("k,kmer_size","k for local assembly",cxxopts::value<uint8_t>(k))
                ("c,min_cvg","minimum coverage for dbg kmers",cxxopts::value<int>(min_cvg))
                ("use_tag_cvg","use tag coverage for dbg kmer",cxxopts::value<bool>(use_tag_cvg));

        options.add_options("Development")
                ("dev_benchmark","file with a benchmark definition and problem list",cxxopts::value<std::string>(benchmark_filename))
                ("dev_problem_analysis","run a LOT of problem analysis... really for debugging, and dump with this prefix",cxxopts::value<std::string>(problem_analysis));

        auto result(options.parse(argc, argv));

        if (result.count("help"))
        {
            std::cout << options.help({"","Assembly","Development"}) << std::endl;
            exit(0);
        }

        if ((result.count("p")!=1 and result.count("dev_benchmarck")) or result.count("o")!=1) {
            throw cxxopts::OptionException(" please specify input problem and output prefix");
        }



    } catch (const cxxopts::OptionException& e)
    {
        std::cout << "Error parsing options: " << e.what() << std::endl << std::endl
                <<"Use option --help to check command line arguments." << std::endl;
        exit(1);
    }

    std::cout<<std::endl;


    //======= BENCHMARK MODE
    std::vector<std::pair<std::string,std::string>> filenames;
    if (!benchmark_filename.empty()){
        std::ifstream benchfile(benchmark_filename);
        while (!benchfile.eof()){
            std::string group,pfile;
            benchfile>>group;
            benchfile>>pfile;
            if (group.empty() or pfile.empty()) break;
            filenames.emplace_back(group, pfile);
        }

        //read all filenames, run each one with all alternatives and record their result stats: n50, anchors fully resolved, patches created and time spent.
        sglib::OutputLog()<<"=== RUNNING BENCHMARK ON "<<filenames.size()<<" PROBLEMS FROM FILE ==="<<std::endl;
        std::ofstream stats_file(output_prefix+"_benchmark_stats.csv");
        std::vector<std::string> pass_names={"Unclean","Clean"};
        for (auto li=0;li<filenames.size();++li) {
            std::vector<std::pair<std::string, std::string>> metrics[pass_names.size()];
            for (auto pass = 0; pass < pass_names.size(); ++pass) {
                //auto &l = lines[li];
                WorkSpace ws;
                LocalHaplotypeAssembler lha(ws);
                lha.init_from_full_file(filenames[li].second);
                if (pass==0) lha.write_anchors(output_prefix+"_p"+std::to_string(li)+"_anchors.fasta");
                std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
                sglib::OutputLog() << "Running \""+pass_names[pass]+"\" on "+filenames[li].second<<std::endl;
                //=== Definition of the passes!!!!
                if (0 == pass) {
                    lha.assemble(63, 5, false, false);
                } else if (1 == pass) {
                    lha.assemble(63, 7, false, false);
                    lha.create_63mer_index();
                    lha.path_linked_reads_informative_singles();
                    lha.expand_canonical_repeats();
                    lha.assembly.join_all_unitigs();
                    lha.create_63mer_index();
                    lha.path_linked_reads_informative_singles();
                    lha.expand_canonical_repeats();
                    lha.assembly.join_all_unitigs();
                    //lha.expand_canonical_repeats_direct(140);
//                    lha.path_linked_reads_informative_singles();
//                    lha.path_paired_reads_informative_singles();
//                    uint64_t pos=0,neg=0;
//                    uint64_t ppos=0,pneg=0;
//                    for (auto kp:lha.assembly.kmer_to_graphposition) {
//                        if (kp.second.node>0) pos++;
//                        else neg++;
//                        if (kp.second.pos>0) ppos++;
//                        else pneg++;
//                    }
//                    std::cout<<"Index has "<<pos<<" positive nodes and "<<neg<<" negatives"<<std::endl;
//                    std::cout<<"Index has "<<ppos<<" positive positions and "<<pneg<<" negatives"<<std::endl;
                    std::ofstream pathsf(output_prefix+"_p"+std::to_string(li)+"_paths.txt");
                    for (auto &p:lha.linkedread_paths){
                        for (auto n:p) pathsf<<" "<<n;
                        pathsf<<std::endl;
                    }


                }
                //=== ENDS
                std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
                metrics[pass] = lha.compute_metrics();
                metrics[pass].emplace_back("Runtime", std::to_string(
                        1.0 * std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count() / 1000));
                lha.write_gfa(output_prefix+"_p"+std::to_string(li)+"_"+pass_names[pass]+".gfa");
            }
            if (li==0) {
                stats_file<<"Filename";
                for (auto pass = 0; pass < pass_names.size(); ++pass)
                    for (auto m:metrics[pass]) stats_file<<","<<pass_names[pass]<<"_"<<m.first;
                stats_file<<std::endl;
            }
            stats_file<<filenames[li].second;
            for (auto pass = 0; pass < pass_names.size(); ++pass)
                for (auto m:metrics[pass]) stats_file<<","<<m.second;
            stats_file<<std::endl;

        }
        //output all the stats in a graph format (maybe gnuplot or something like that? /3ds?
        exit(0);
    }

    //======= WORKSPACE LOAD AND CHECKS ======
    WorkSpace ws;
    LocalHaplotypeAssembler lha(ws);

    if (!workspace_file.empty()) {
        sglib::OutputLog() << "Loading Workspace..." << std::endl;
        ws.load_from_disk(workspace_file);
        if (!ws.sg.is_sane()) {
            sglib::OutputLog() << "ERROR: sg.is_sane() = false" << std::endl;
            return 1;
        }

        //TODO: other checks? reads mapped to valid nodes and such?

        sglib::OutputLog() << "Loading Workspace DONE" << std::endl;
        lha.init_from_file(problem_file);
    }
    else {
        lha.init_from_full_file(problem_file);
    }
    if (!problem_analysis.empty()) {
        lha.problem_analysis(problem_analysis);
        exit(0);
    }
    //TODO: try different coverage cutoffs and keep the patches from all of them
    lha.assemble(63, min_cvg, use_tag_cvg, false);
    lha.create_63mer_index();
    lha.path_linked_reads_informative_singles();
    lha.expand_canonical_repeats();
    lha.assembly.join_all_unitigs();
    lha.create_63mer_index();
    lha.path_linked_reads_informative_singles();
    lha.expand_canonical_repeats();
    lha.assembly.join_all_unitigs();
    lha.write_anchors(output_prefix+"_anchors.fasta");
    lha.write_gfa(output_prefix+"_final.gfa");
    lha.construct_patches();
    lha.write_patches(output_prefix+"_patches.fasta");
    lha.construct_patched_backbone();
    lha.write_patched_backbone(output_prefix+"_patchedbackbone.fasta");
    lha.construct_patched_backbone(false);
    lha.write_patched_backbone(output_prefix+"_patchedbackbone_multi.fasta");

    return 0;
}

