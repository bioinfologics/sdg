#include <vector>
#include <iostream>
#include <fstream>
#include <sdglib/processors/KmerCompressionIndex.hpp>
#include <sdglib/workspace/WorkSpace.hpp>
#include "cxxopts.hpp"

enum class WorkspaceFunctions{
    NONE, //This is needed because the default of a map is 0
    MAKE,
    STATUS,
    DUMP,
    ADD_DS
};
struct WorkspaceFunctionMap : public std::map<std::string, WorkspaceFunctions>
{
    WorkspaceFunctionMap()
    {
        operator[]("make") =  WorkspaceFunctions::MAKE;
        operator[]("status") = WorkspaceFunctions::STATUS;
        operator[]("dump") = WorkspaceFunctions::DUMP;
        operator[]("add") = WorkspaceFunctions::ADD_DS;
    };
};

void make_workspace(int argc, char** argv){
    std::vector<std::string> lr_datastores, pr_datastores, Lr_datastores;
    std::string output;
    std::string gfa_filename;
    std::string kci_filename;
    std::string ws_filename;
    try {
        cxxopts::Options options("sdg-workspace make", "SDG make workspace");

        options.add_options()
                ("help", "Print help")
                ("g,gfa", "input gfa file", cxxopts::value(gfa_filename))
                ("w,workspace", "input workspace file", cxxopts::value(ws_filename))
                ("p,paired_reads", "paired reads datastore", cxxopts::value(pr_datastores))
                ("l,linked_reads", "linked reads datastore", cxxopts::value(lr_datastores))
                ("L,long_reads", "long reads datastore", cxxopts::value(Lr_datastores))
                ("o,output", "output file", cxxopts::value(output));
        auto newargc = argc - 1;
        auto newargv = &argv[1];
        auto result = options.parse(newargc, newargv);
        if (result.count("help")) {
            std::cout << options.help({""}) << std::endl;
            exit(0);
        }

        if (!output.empty()) {
            throw cxxopts::OptionException(" please specify an output prefix");
        }
        if ( (gfa_filename.empty() and ws_filename.empty()) or (!gfa_filename.empty() and !ws_filename.empty() ) ) {
            throw cxxopts::OptionException(" please specify either a gfa or workspace input file");
        }

    } catch (const cxxopts::OptionException &e) {
        std::cout << "Error parsing options: " << e.what() << std::endl << std::endl
                  << "Use option --help to check command line arguments." << std::endl;
        exit(1);
    }

    //===== LOAD GRAPH =====
    WorkSpace output_ws;
    OperationJournal op;
    if (!ws_filename.empty()) {
        WorkSpace w;
        w.load_from_disk(ws_filename);
        output_ws.linked_reads_datastores = w.linked_reads_datastores;
        output_ws.paired_reads_datastores = w.paired_reads_datastores;
        output_ws.long_reads_datastores = w.long_reads_datastores;
        output_ws.operation_journals = w.operation_journals;
        output_ws.kmer_counts_datastores = w.kmer_counts_datastores;
        output_ws.sdg = w.sdg;
        output_ws.distance_graphs = w.distance_graphs;
        op = output_ws.add_operation("Copy", std::string("sdg-workspace make ") + GIT_ORIGIN_URL + " -> " + GIT_BRANCH + " " +
                                        GIT_COMMIT_HASH, std::string("Copied with sdg-workspace make from")+ws_filename);
    }
    if (!gfa_filename.empty()) {
        op = output_ws.add_operation("Creation",
                                           std::string("sdg-workspace make ") + GIT_ORIGIN_URL + " -> " + GIT_BRANCH +
                                           " " +
                                           GIT_COMMIT_HASH, "Created with sdg-workspace make");
        output_ws.sdg.load_from_gfa(gfa_filename);
        op.addEntry("GFA imported from " + gfa_filename + " (" + std::to_string(output_ws.sdg.nodes.size() - 1) +
                     " nodes)");
    }

    for (auto prds:pr_datastores) {
        //create and load the datastore, and the mapper!
        output_ws.paired_reads_datastores.emplace_back(output_ws, prds);
        op.addEntry("PairedReadDatastore imported from " + prds + " (" +
                     std::to_string(output_ws.paired_reads_datastores.back().size()) + " reads)");
    }

    for (auto lrds:lr_datastores) {
        //create and load the datastore, and the mapper!
        output_ws.linked_reads_datastores.emplace_back(output_ws, lrds);
        op.addEntry("LinkedReadDatastore imported from " + lrds + " (" +
                     std::to_string(output_ws.linked_reads_datastores.back().size()) + " reads)");
    }

    for (auto Lrds:Lr_datastores) {
        //create and load the datastore, and the mapper!
        output_ws.long_reads_datastores.emplace_back(output_ws, Lrds);
        op.addEntry("LongReadDatastore imported from " + Lrds + " (" +
                     std::to_string(output_ws.long_reads_datastores.back().size()) + " reads)");
    }

    output_ws.dump_to_disk(output + ".bsgws");
}

void status_workspace(int argc, char **argv){
    std::string filename;
    try {

        cxxopts::Options options("sdg-workspace status", "SDG workspace status");

        options.add_options()
                ("help", "Print help")
                ("w,workspace", "workspace filename", cxxopts::value<std::string>(filename));

        auto newargc=argc-1;
        auto newargv=&argv[1];
        auto result=options.parse(newargc,newargv);
        if (result.count("help")) {
            std::cout << options.help({""}) << std::endl;
            exit(0);
        }

        if (result.count("workspace")==0) {
            throw cxxopts::OptionException(" please specify workspace file");
        }


    } catch (const cxxopts::OptionException &e) {
        std::cout << "Error parsing options: " << e.what() << std::endl << std::endl
                  << "Use option --help to check command line arguments." << std::endl;
        exit(1);
    }
    WorkSpace w;
    w.load_from_disk(filename);
    w.status();
    std::cout<<std::endl<<"---=== Workspace current status ===---"<<std::endl;
}

void dump_workspace(int argc, char **argv){
    std::string filename,gfafilename;
    try {

        cxxopts::Options options("sdg-workspace dump", "SDG workspace dump");

        options.add_options()
                ("help", "Print help")
                ("w,workspace", "workspace filename", cxxopts::value<std::string>(filename))
                ("g,gfa", "gfa output prefix", cxxopts::value<std::string>(gfafilename));

        auto newargc=argc-1;
        auto newargv=&argv[1];
        auto result=options.parse(newargc,newargv);
        if (result.count("help")) {
            std::cout << options.help({""}) << std::endl;
            exit(0);
        }

        if (result.count("workspace")==0) {
            throw cxxopts::OptionException(" please specify kmer spectra file");
        }


    } catch (const cxxopts::OptionException &e) {
        std::cout << "Error parsing options: " << e.what() << std::endl << std::endl
                  << "Use option --help to check command line arguments." << std::endl;
        exit(1);
    }
    std::cout << "Loading workspace " << std::endl;
    WorkSpace w;
    w.load_from_disk(filename);
    if (!w.sdg.is_sane()) {
        sdglib::OutputLog()<<"ERROR: sdg.is_sane() = false"<<std::endl;
        //return 1;
    }
    std::cout << "Done ... " << std::endl;
    std::cout << "Dumping gfa workspace " << std::endl;
    if (not gfafilename.empty()){
        w.sdg.write_to_gfa1(gfafilename + ".gfa");
    }
    std::cout << "Done... " << std::endl;
}

int main(int argc, char * argv[]) {
    std::cout << "sdg-workspace"<<std::endl<<std::endl;
    std::cout << "Git origin: " << GIT_ORIGIN_URL << " -> "  << GIT_BRANCH << std::endl;
    std::cout << "Git commit: " << GIT_COMMIT_HASH << std::endl<<std::endl;
    std::cout << "Executed command:"<<std::endl;

    WorkspaceFunctionMap functionMap;

    for (auto i=0;i<argc;i++) std::cout<<argv[i]<<" ";
    std::cout<<std::endl<<std::endl;
    if (argc <2){
        std::cout<<"Please specify one of the following functions: ";
        for (auto &p: functionMap) {std::cout << p.first << ", ";}
        std::cout << "to select the execution mode." << std::endl;
        exit(1);
    }

    auto function = functionMap.find(argv[1])->second;

    switch(function) {
        case WorkspaceFunctions::MAKE:
            make_workspace(argc,argv);
            break;
        case WorkspaceFunctions::STATUS:
            status_workspace(argc,argv);
            break;
        case WorkspaceFunctions::DUMP:
            dump_workspace(argc,argv);
            break;
        default:
            std::cout << "Invalid option specified." << std::endl;
            std::cout<<"Please specify one of the following functions: ";
            for (auto &p: functionMap) {std::cout << p.first << ", ";}
            std::cout << "to select the execution mode." << std::endl;
            exit(1);
            break;
    }
    exit(0);
}
