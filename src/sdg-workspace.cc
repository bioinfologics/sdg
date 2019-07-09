#include <vector>
#include <iostream>
#include <fstream>
#include <sdglib/processors/KmerCompressionIndex.hpp>
#include <sdglib/workspace/WorkSpace.hpp>
#include "cxxopts.hpp"

enum class WorkspaceFunctions{
    NONE, //This is needed because the default of a map is 0
    MAKE,
    LOG,
    DUMP,
    ADD_DS
};
struct WorkspaceFunctionMap : public std::map<std::string, WorkspaceFunctions>
{
    WorkspaceFunctionMap()
    {
        operator[]("make") =  WorkspaceFunctions::MAKE;
        operator[]("log") = WorkspaceFunctions::LOG;
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

void log_workspace(int argc, char **argv){
    std::string filename;
    try {

        cxxopts::Options options("sdg-workspace log", "SDG workspace log");

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
    //graph
    w.sdg.print_status();

    //PR datastores and mappings
    sdglib::OutputLog()<<"Workspace contains "<< w.paired_reads_datastores.size() << " paired reads datastores" <<std::endl;
    for (auto di=0;di<w.paired_reads_datastores.size();++di){
        w.paired_reads_datastores[di].print_status();

    }
    //10x datastores and mappings
    sdglib::OutputLog()<<"Workspace contains "<< w.linked_reads_datastores.size() << " linked reads datastores" <<std::endl;
    for (auto di=0;di<w.linked_reads_datastores.size();++di){
        w.linked_reads_datastores[di].print_status();
    }
    //LR datastores and mappings
    sdglib::OutputLog()<<"Workspace contains "<< w.long_reads_datastores.size() << " long reads datastores" <<std::endl;
    for (auto di=0;di<w.long_reads_datastores.size();++di){
        w.long_reads_datastores[di].print_status();
    }
}

void dump_workspace(int argc, char **argv){
    std::string filename,gfafilename,nodeinfofilename,seqfilename;
    float minKCI=.5, maxKCI=1.25;
    size_t min_size=2000,max_size=100000000;
    try {

        cxxopts::Options options("sdg-workspace dump", "SDG workspace dump");

        options.add_options()
                ("help", "Print help")
                ("w,workspace", "workspace filename", cxxopts::value<std::string>(filename))
                ("g,gfa", "gfa output prefix", cxxopts::value<std::string>(gfafilename))
                ("n,node_info", "node info prefix",cxxopts::value<std::string>(nodeinfofilename))
                ("s,node_sequences", "selected sequences prefix", cxxopts::value<std::string>(seqfilename))
                ("min_kci", "selected sequences min KCI", cxxopts::value<float>(minKCI))
                ("max_kci", "selected sequences max KCI", cxxopts::value<float>(maxKCI))
                ("min_size", "selected sequences min size", cxxopts::value<size_t>(min_size))
                ("max_size", "selected sequences max size", cxxopts::value<size_t>(max_size));

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

void add_datastores(int argc, char **argv) {
    std::string output;
    std::string base_filename;
    std::vector<std::string> add_filenames;
    try {
        cxxopts::Options options("sdg-workspace add", "SDG workspace add, adds workspaces to a base workspace in the order they are passed in to an output workspace.");

        options.add_options()
                ("help", "Print help")
                ("w,workspace_base", "Base workspace, the graph and KCI are taken from this workspace", cxxopts::value<std::string>(base_filename))
                ("a,workspace_add", "Workspaces to add, can be defined multiple times (they will be added in the order specified)", cxxopts::value<std::vector<std::string>>(add_filenames))
                ("o,output", "Output workspace", cxxopts::value<std::string>(output));
        auto newargc=argc-1;
        auto newargv=&argv[1];
        auto result=options.parse(newargc,newargv);
        if (result.count("help")) {
            std::cout << options.help({""}) << std::endl;
            exit(0);
        }

        if (output=="" or base_filename=="" or add_filenames.empty()) {
            throw cxxopts::OptionException(" please specify a base workspace, a merge workspace and output prefix");
        }

    } catch (const cxxopts::OptionException &e) {
        std::cout << "Error parsing options: " << e.what() << std::endl << std::endl
                  << "Use option --help to check command line arguments." << std::endl;
        exit(1);
    }

    WorkSpace base,out;

    base.load_from_disk(base_filename);
    out.sdg = base.sdg;

    /* Copy BASE datastores and mappers */
    for (int i = 0; i < base.paired_reads_datastores.size(); ++i){
        out.paired_reads_datastores.emplace_back(out, base.paired_reads_datastores[i]);
    }
    for (int i = 0; i < base.linked_reads_datastores.size(); ++i){
        out.linked_reads_datastores.emplace_back(out, base.linked_reads_datastores[i]);
    }
    for (int i = 0; i < base.long_reads_datastores.size(); ++i){
        out.long_reads_datastores.emplace_back(out, base.long_reads_datastores[i]);
    }

    /* Copy for each ADD their datastores and mappers */
    for (const auto &fn: add_filenames) {
        WorkSpace add;
        add.load_from_disk(fn);
        for (int i = 0; i < add.paired_reads_datastores.size(); ++i){
            out.paired_reads_datastores.emplace_back(out, add.paired_reads_datastores[i]);
        }
        for (int i = 0; i < add.linked_reads_datastores.size(); ++i){
            out.linked_reads_datastores.emplace_back(out, add.linked_reads_datastores[i]);
        }
        for (int i = 0; i < add.long_reads_datastores.size(); ++i){
            out.long_reads_datastores.emplace_back(out, add.long_reads_datastores[i]);
        }
    }
    out.dump_to_disk(output+".bsgws");
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
        case WorkspaceFunctions::LOG:
            log_workspace(argc,argv);
            break;
        case WorkspaceFunctions::DUMP:
            dump_workspace(argc,argv);
            break;
        case WorkspaceFunctions::ADD_DS:
            add_datastores(argc,argv);
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
