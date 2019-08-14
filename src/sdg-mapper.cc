#include <iostream>
#include <fstream>
#include <sdglib/workspace/WorkSpace.hpp>
#include "sdglib/utilities/OutputLog.hpp"
#include "cxxopts.hpp"

int main(int argc, char * argv[]) {
    std::cout << "Welcome to sdg-mapper"<<std::endl<<std::endl;
    std::cout << "Git origin: " << GIT_ORIGIN_URL << " -> "  << GIT_BRANCH << std::endl;
    std::cout << "Git commit: " << GIT_COMMIT_HASH << std::endl<<std::endl;
    std::cout << "Executed command:"<<std::endl;
    bool sat_kmer_index=false;
    unsigned int long_reads_k=15;
    bool use63mers=false,best_nodes=false;
    bool skip_paired=false,skip_linked=false,skip_long=false;
    unsigned int max_filter=200;
    for (auto i=0;i<argc;i++) std::cout<<argv[i]<<" ";
    std::cout<<std::endl<<std::endl;

    std::string workspace_file,output_prefix;
    sdglib::OutputLogLevel=sdglib::LogLevels::DEBUG;
    try
    {
        cxxopts::Options options("sdg-mapper", "reads-to-graph mapper for sdg worskpaces");

        options.add_options()
                ("help", "Print help")
                ("w,workspace", "input workspace", cxxopts::value<std::string>(workspace_file))
                ("o,output", "output file prefix", cxxopts::value<std::string>(output_prefix))
                ("k","long read indexing/mapping kmer size", cxxopts::value(long_reads_k)->default_value("15"))
                ("s,use_sat-index", "Use saturated small-k index", cxxopts::value(sat_kmer_index)->default_value("false")->implicit_value("true"))
                ("b,best_nodes_only", "Map long reads to best nodes only", cxxopts::value(best_nodes))
                ("m,max_kmer_repeat", "maximum number of times a kmer appears (LongReadMapper)", cxxopts::value(max_filter)->default_value("200"))
                ("use_63-mers", "mapping based on 63-mers", cxxopts::value<bool>(use63mers))
                ("skip_paired", "skip paired read mapping", cxxopts::value<bool>(skip_paired))
                ("skip_linked", "skip paired read mapping", cxxopts::value<bool>(skip_linked))
                ("skip_long", "skip paired read mapping", cxxopts::value<bool>(skip_long));



        auto result(options.parse(argc, argv));

        if (result.count("help"))
        {
            std::cout << options.help({""}) << std::endl;
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
    sdglib::OutputLog()<<"Loading Workspace..."<<std::endl;
    ws.load_from_disk(workspace_file);
    auto op = ws.add_operation("Mapping reads", std::string("sdg-mapper") + GIT_ORIGIN_URL + " -> " + GIT_BRANCH +
                                      " " +
                                      GIT_COMMIT_HASH, "Mapping reads");
    op.addEntry("sdg-mapper run started");
    sdglib::OutputLog()<<"Loading Workspace DONE"<<std::endl;
    sdglib::OutputLog()<<"Mapping reads..."<<std::endl;
    if (not skip_paired) {
        auto pri = 0;
        for (auto &ds:ws.paired_reads_datastores) {
            sdglib::OutputLog() << "Mapping reads from paired library..." << std::endl;
            if (!use63mers) ds.mapper.remap_all_reads();
            else ds.mapper.remap_all_reads63();
            ds.mapper.print_status();
            sdglib::OutputLog() << "Computing size distribution..." << std::endl;
            auto sdist = ds.mapper.size_distribution();
            std::ofstream df("prdist_" + std::to_string(pri++) + ".csv");
            for (auto i = 0; i < sdist.size(); i += 10) {
                uint64_t t = 0;
                for (auto j = i; j < i + 10; ++j) t += sdist[j];
                if (t > 0) df << i << ", " << t << std::endl;
            }
            op.addEntry("reads from " + ds.filename + " re-mapped to current graph");
            sdglib::OutputLog() << "Mapping reads from paired library DONE." << std::endl;
        }
    }
    if (not skip_linked) {
        for (auto &ds:ws.linked_reads_datastores) {
            sdglib::OutputLog() << "Mapping reads from linked library..." << std::endl;
            if (!use63mers) ds.mapper.remap_all_reads();
            else ds.mapper.remap_all_reads63();
            op.addEntry("reads from " + ds.filename + " re-mapped to current graph");
            sdglib::OutputLog() << "Mapping reads from linked library DONE." << std::endl;
        }
    }
    if (not skip_long) {
        for (auto &ds: ws.long_reads_datastores) {
            sdglib::OutputLog() << "Mapping reads from long reads library..." << std::endl;
            ds.mapper.sat_kmer_index = sat_kmer_index;
            ds.mapper.k = long_reads_k;
            ds.mapper.max_index_freq = max_filter;
            if (not best_nodes) ds.mapper.map_reads();
            else ds.mapper.map_reads_to_best_nodes();
            op.addEntry("reads from " + ds.filename + " re-mapped to current graph");
            sdglib::OutputLog() << "Mapping reads from long reads library DONE." << std::endl;
        }
    }
    op.addEntry("sdg-mapper run finished");
    ws.dump_to_disk(output_prefix+".sdgws");
    sdglib::OutputLog()<<"Mapping reads DONE."<<std::endl;
    return 0;
}

