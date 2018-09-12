#include <iostream>
#include <fstream>
#include <sglib/datastores/PairedReadsDatastore.hpp>
#include <sglib/processors/GraphMaker.hpp>
#include <sglib/workspace/WorkSpace.hpp>
#include "cxxopts.hpp"


int main(int argc, char * argv[]) {
    std::cout << "Welcome to bsg-dbg"<<std::endl<<std::endl;
    std::cout << "Git origin: " << GIT_ORIGIN_URL << " -> "  << GIT_BRANCH << std::endl;
    std::cout << "Git commit: " << GIT_COMMIT_HASH << std::endl<<std::endl;
    std::cout << "Executed command:"<<std::endl;
    for (auto i=0;i<argc;i++) std::cout<<argv[i]<<" ";
    std::cout<<std::endl<<std::endl;

    std::string pr_file,lr_file,output_prefix;
    int min_coverage=5,k=63;
    sglib::OutputLogLevel=sglib::LogLevels::INFO;
    try
    {
        cxxopts::Options options("bsg-dbg", "create a DBG graph from a short-read datastore, and populate KCI");

        options.add_options()
                ("help", "Print help")
                ("p,paired_datastore", "input paired read datastore", cxxopts::value<std::string>(pr_file))
                //("l,linked_datastore", "input linked read datastore", cxxopts::value<std::string>(lr_filename))
                ("o,output", "output file prefix", cxxopts::value<std::string>(output_prefix));
        options.add_options("Heuristics")
                ("k","k value for DBG construction (default: 63)",cxxopts::value<int>(k))
                ("c,min_coverage","minimum coverage for a kmer to include in DBG",cxxopts::value<int>(min_coverage))
                ;




        auto result(options.parse(argc, argv));

        if (result.count("help"))
        {
            std::cout << options.help({"","Heuristics"}) << std::endl;
            exit(0);
        }

        if (result.count("p")!=1 or result.count("o")!=1) {
            throw cxxopts::OptionException(" please specify input datastore and output prefix");
        }



    } catch (const cxxopts::OptionException& e)
    {
        std::cout << "Error parsing options: " << e.what() << std::endl << std::endl
                  <<"Use option --help to check command line arguments." << std::endl;
        exit(1);
    }

    std::cout<<std::endl;
    WorkSpace ws;
    ws.add_log_entry("Workspace created with bsg-dbg");
    ws.add_log_entry("Origin datastore: "+pr_file);
    ws.paired_read_datastores.emplace_back(pr_file);
    ws.paired_read_mappers.emplace_back(ws.sg,ws.paired_read_datastores.back(), ws.uniqueKmerIndex, ws.unique63merIndex);
    GraphMaker gm(ws.sg);
    //counting kmers...
    sglib::OutputLog()<<"Creating "<<k<<"-mer set from datastore..."<<std::endl;
    auto kmers=ws.paired_read_datastores.back().get_all_kmers128(k,min_coverage);
    sglib::OutputLog()<<"DONE! "<<kmers.size()<<" "<<k<<"-mers with coverage >="<<min_coverage<<std::endl;
    sglib::OutputLog()<<"Creating DBG..."<<std::endl;
    gm.new_graph_from_kmerset_trivial128(kmers,k);
    sglib::OutputLog()<<"DONE! "<<ws.sg.count_active_nodes()<<" nodes in graph"<<std::endl;

    gm.tip_clipping(200);

    sglib::OutputLog()<<"DONE! "<<ws.sg.count_active_nodes()<<" nodes in graph"<<std::endl;
    ws.kci.index_graph();
    ws.kci.start_new_count();
    ws.kci.add_counts_from_datastore(ws.paired_read_datastores.back());

    ws.sg.write_to_gfa(output_prefix+"_DBG.gfa");
    ws.dump_to_disk(output_prefix+".bsgws");
}

