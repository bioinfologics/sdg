#include <iostream>
#include <fstream>
#include <sdglib/workspace/WorkSpace.hpp>
#include "sdglib/logger/OutputLog.hpp"
#include "cxxopts.hpp"

int main(int argc, char * argv[]) {
    std::cout << "Welcome to bsg-mapper"<<std::endl<<std::endl;
    std::cout << "Git origin: " << GIT_ORIGIN_URL << " -> "  << GIT_BRANCH << std::endl;
    std::cout << "Git commit: " << GIT_COMMIT_HASH << std::endl<<std::endl;
    std::cout << "Executed command:"<<std::endl;
    bool use63mers=false;
    unsigned int max_filter=200;
    for (auto i=0;i<argc;i++) std::cout<<argv[i]<<" ";
    std::cout<<std::endl<<std::endl;

    std::string workspace_file,output_prefix;
    sdglib::OutputLogLevel=sdglib::LogLevels::DEBUG;
    try
    {
        cxxopts::Options options("bsg-mapper", "reads-to-graph mapper for bsg worskpaces");

        options.add_options()
                ("help", "Print help")
                ("w,workspace", "input workspace", cxxopts::value<std::string>(workspace_file))
                ("o,output", "output file prefix", cxxopts::value<std::string>(output_prefix))
                ("m,max_kmer_repeat", "maximum number of times a kmer appears (LongReadMapper)", cxxopts::value(max_filter)->default_value("200"))
                ("use_63-mers", "mapping based on 63-mers", cxxopts::value<bool>(use63mers));



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
    ws.add_log_entry("bsg-mapper run started");
    sdglib::OutputLog()<<"Loading Workspace DONE"<<std::endl;
    sdglib::OutputLog()<<"Mapping reads..."<<std::endl;
    auto pri=0;
    if (!ws.paired_read_datastores.empty() or !ws.linked_read_datastores.empty()) {
        if (!use63mers) ws.create_index();
        else ws.create_63mer_index();
    }
    for (auto &m:ws.paired_read_mappers) {
        sdglib::OutputLog()<<"Mapping reads from paired library..."<<std::endl;
        if (!use63mers) m.remap_all_reads();
        else m.remap_all_reads63();
        m.print_status();
        sdglib::OutputLog()<<"Computing size distribution..."<<std::endl;
        auto sdist=m.size_distribution();
        std::ofstream df("prdist_"+std::to_string(pri++)+".csv");
        for (auto i=0;i<sdist.size();i+=10){
            uint64_t t=0;
            for (auto j=i;j<i+10;++j) t+=sdist[j];
            if (t>0) df<<i<<", "<<t<<std::endl;
        }
        ws.add_log_entry("reads from "+m.datastore.filename+" re-mapped to current graph");
        sdglib::OutputLog()<<"Mapping reads from paired library DONE."<<std::endl;
    }
    for (auto &m:ws.linked_read_mappers) {
        sdglib::OutputLog()<<"Mapping reads from linked library..."<<std::endl;
        if (!use63mers) m.remap_all_reads();
        else m.remap_all_reads63();
        ws.add_log_entry("reads from "+m.datastore.filename+" re-mapped to current graph");
        sdglib::OutputLog()<<"Mapping reads from linked library DONE."<<std::endl;
    }
    for (auto &m: ws.long_read_mappers) {
        sdglib::OutputLog()<<"Mapping reads from long reads library..."<<std::endl;
        m.map_reads(max_filter);
        ws.add_log_entry("reads from "+m.datastore.filename+" re-mapped to current graph");
        sdglib::OutputLog()<<"Mapping reads from long reads library DONE."<<std::endl;
    }
    ws.add_log_entry("bsg-mapper run finished");
    ws.dump_to_disk(output_prefix+".bsgws");
    sdglib::OutputLog()<<"Mapping reads DONE."<<std::endl;
    return 0;
}

