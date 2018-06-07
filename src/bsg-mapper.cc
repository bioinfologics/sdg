#include <iostream>
#include <fstream>
#include <sglib/WorkSpace.hpp>
#include <sglib/processors/Untangler.hpp>
#include "sglib/logger/OutputLog.h"
#include "cxxopts.hpp"

int main(int argc, char * argv[]) {
    std::cout << "Welcome to bsg-mapper"<<std::endl<<std::endl;
    std::cout << "Git origin: " << GIT_ORIGIN_URL << " -> "  << GIT_BRANCH << std::endl;
    std::cout << "Git commit: " << GIT_COMMIT_HASH << std::endl<<std::endl;
    std::cout << "Executed command:"<<std::endl;
    for (auto i=0;i<argc;i++) std::cout<<argv[i]<<" ";
    std::cout<<std::endl<<std::endl;

    std::string workspace_file,output_prefix;
    sglib::OutputLogLevel=sglib::LogLevels::DEBUG;
    try
    {
        cxxopts::Options options("bsg-mapper", "reads-to-graph mapper for bsg worskpaces");

        options.add_options()
                ("help", "Print help")
                ("w,workspace", "input workspace", cxxopts::value<std::string>(workspace_file))
                ("o,output", "output file prefix", cxxopts::value<std::string>(output_prefix));



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
    sglib::OutputLog()<<"Loading Workspace..."<<std::endl;
    ws.load_from_disk(workspace_file);
    ws.add_log_entry("bsg-mapper run started");
    sglib::OutputLog()<<"Loading Workspace DONE"<<std::endl;
    sglib::OutputLog()<<"Mapping reads..."<<std::endl;
    auto pri=0;
    ws.getGraph().create_index();
    for (auto &m:ws.getPairedReadMappers()) {
        sglib::OutputLog()<<"Mapping reads from paired library..."<<std::endl;
        m.map_reads();
        m.print_stats();
        sglib::OutputLog()<<"Computing size distribution..."<<std::endl;
        auto sdist=m.size_distribution();
        std::ofstream df("prdist_"+std::to_string(pri++)+".csv");
        for (auto i=0;i<sdist.size();i+=10){
            uint64_t t=0;
            for (auto j=i;j<i+10;++j) t+=sdist[j];
            if (t>0) df<<i<<", "<<t<<std::endl;
        }
        ws.add_log_entry("reads from "+m.datastore.filename+" re-mapped to current graph");
        sglib::OutputLog()<<"Mapping reads from paired library DONE."<<std::endl;
    }
    for (auto &m:ws.getPairedReadMappers()) {
        sglib::OutputLog()<<"Mapping reads from linked library..."<<std::endl;
        m.map_reads();
        ws.add_log_entry("reads from "+m.datastore.filename+" re-mapped to current graph");
        sglib::OutputLog()<<"Mapping reads from linked library DONE."<<std::endl;
    }
    ws.getPathsDatastore().clear();
    ws.add_log_entry("path_datastores cleared");
    ws.add_log_entry("bsg-mapper run finished");
    ws.dump_to_disk(output_prefix+".bsgws");
    sglib::OutputLog()<<"Mapping reads DONE."<<std::endl;
    return 0;
}

