#include <vector>
#include <iostream>
#include <fstream>
#include <sdglib/workspace/WorkSpace.hpp>
#include "cxxopts.hpp"

enum class WorkspaceFunctions{
    NONE, //This is needed because the default of a map is 0
    MAKE,
    STATUS,
    DUMP,
    ADD_DS,
    ADD_COUNTS
};
struct WorkspaceFunctionMap : public std::map<std::string, WorkspaceFunctions>
{
    WorkspaceFunctionMap()
    {
        operator[]("make") =  WorkspaceFunctions::MAKE;
        operator[]("status") = WorkspaceFunctions::STATUS;
        operator[]("dump") = WorkspaceFunctions::DUMP;
        operator[]("add_ds") = WorkspaceFunctions::ADD_DS;
        operator[]("add_counts") = WorkspaceFunctions::ADD_COUNTS;
    };
};


void make_workspace(int argc, char** argv){
    std::string toolname("sdg-workspace make: ");
    std::string git_version(std::string(GIT_ORIGIN_URL) + " -> " + GIT_BRANCH + " " +
                            GIT_COMMIT_HASH);
    std::vector<std::string> lr_datastores, pr_datastores, Lr_datastores;
    std::string gfa_filename;
    std::string kci_filename;
    std::string ws_filename;
    std::string output;
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

        if (output.empty()) {
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
    if (!ws_filename.empty()) {
        WorkSpace w;
        w.load_from_disk(ws_filename);
        output_ws.linked_reads_datastores = w.linked_reads_datastores;
        output_ws.paired_reads_datastores = w.paired_reads_datastores;
        output_ws.long_reads_datastores = w.long_reads_datastores;
        output_ws.journal = w.journal;
        output_ws.kmer_counts = w.kmer_counts;
        output_ws.sdg = w.sdg;
        output_ws.distance_graphs = w.distance_graphs;
        auto op = output_ws.add_operation("Copy", toolname + git_version, std::string("Copied with sdg-workspace make from ")+ws_filename);
    }
    if (!gfa_filename.empty()) {
        auto op = output_ws.add_operation("Create", toolname + git_version, "Created with sdg-workspace make");
        output_ws.sdg.load_from_gfa(gfa_filename);
        op.addEntry("GFA imported from " + gfa_filename + " (" + std::to_string(output_ws.sdg.nodes.size() - 1) +
                     " nodes)");
    }

    for (auto prds:pr_datastores) {
        //create and load the datastore, and the mapper!
        output_ws.paired_reads_datastores.emplace_back(output_ws, prds);
        output_ws.journal.back().addEntry("PairedReadDatastore imported from " + prds + " (" +
                     std::to_string(output_ws.paired_reads_datastores.back().size()) + " reads)");
    }

    for (auto lrds:lr_datastores) {
        //create and load the datastore, and the mapper!
        output_ws.linked_reads_datastores.emplace_back(output_ws, lrds);
        output_ws.journal.back().addEntry("LinkedReadDatastore imported from " + lrds + " (" +
                     std::to_string(output_ws.linked_reads_datastores.back().size()) + " reads)");
    }

    for (auto Lrds:Lr_datastores) {
        //create and load the datastore, and the mapper!
        output_ws.long_reads_datastores.emplace_back(output_ws, Lrds);
        output_ws.journal.back().addEntry("LongReadDatastore imported from " + Lrds + " (" +
                     std::to_string(output_ws.long_reads_datastores.back().size()) + " reads)");
    }

    output_ws.dump_to_disk(output + ".sdgws");
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

void add_datastore(int argc, char **argv) {
    std::string toolname("sdg-workspace add_ds: ");
    std::string git_version(std::string(GIT_ORIGIN_URL) + " -> " + GIT_BRANCH + " " +
                            GIT_COMMIT_HASH);
    std::vector<std::string> lr_datastores, pr_datastores, Lr_datastores;
    std::string ws_filename;
    std::string output;

    try {
        cxxopts::Options options("sdg-workspace add_ds", "SDG add datastore to workspace");

        options.add_options()
                ("help", "Print help")
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

        if (output.empty()) {
            throw cxxopts::OptionException(" please specify an output prefix");
        }
        if ( ws_filename.empty() ) {
            throw cxxopts::OptionException(" please specify a workspace input file");
        }

    } catch (const cxxopts::OptionException &e) {
        std::cout << "Error parsing options: " << e.what() << std::endl << std::endl
                  << "Use option --help to check command line arguments." << std::endl;
        exit(1);
    }
    WorkSpace output_ws(ws_filename);
    output_ws.add_operation("ADD_DS", toolname + git_version, std::string("Added datastores with sdg-workspace add_ds from ")+ws_filename);
    for (auto prds:pr_datastores) {
        output_ws.paired_reads_datastores.emplace_back(output_ws, prds);
        output_ws.journal.back().addEntry("PairedReadDatastore imported from " + prds + " (" +
                    std::to_string(output_ws.paired_reads_datastores.back().size()) + " reads)");
    }

    for (auto lrds:lr_datastores) {
        output_ws.linked_reads_datastores.emplace_back(output_ws, lrds);
        output_ws.journal.back().addEntry("LinkedReadDatastore imported from " + lrds + " (" +
                    std::to_string(output_ws.linked_reads_datastores.back().size()) + " reads)");
    }

    for (auto Lrds:Lr_datastores) {
        output_ws.long_reads_datastores.emplace_back(output_ws, Lrds);
        output_ws.journal.back().addEntry("LongReadDatastore imported from " + Lrds + " (" +
                    std::to_string(output_ws.long_reads_datastores.back().size()) + " reads)");
    }
    output_ws.dump_to_disk(output + ".bsgws");
}

void add_counts(int argc, char **argv) {
    std::string toolname("sdg-workspace add_counts: ");
    std::string git_version(std::string(GIT_ORIGIN_URL) + " -> " + GIT_BRANCH + " " +
                            GIT_COMMIT_HASH);
    std::string ws_filename;
    std::string counts_filename;
    std::string name;
    std::string output;

    try {
        cxxopts::Options options("sdg-workspace add_counts", "SDG add KmerCounts to workspace");

        options.add_options()
                ("help", "Print help")
                ("w,workspace", "input workspace file", cxxopts::value(ws_filename))
                ("n,name", "counts name", cxxopts::value(name))
                ("c,counts_file", "input counts file", cxxopts::value(counts_filename))
                ("o,output", "output file", cxxopts::value(output));
        auto newargc = argc - 1;
        auto newargv = &argv[1];
        auto result = options.parse(newargc, newargv);
        if (result.count("help")) {
            std::cout << options.help({""}) << std::endl;
            exit(0);
        }

        if (output.empty()) {
            throw cxxopts::OptionException(" please specify an output prefix");
        }
        if ( ws_filename.empty() ) {
            throw cxxopts::OptionException(" please specify a workspace input file");
        }

    } catch (const cxxopts::OptionException &e) {
        std::cout << "Error parsing options: " << e.what() << std::endl << std::endl
                  << "Use option --help to check command line arguments." << std::endl;
        exit(1);
    }
    WorkSpace output_ws(ws_filename);
    output_ws.add_operation("ADD_COUNTS", toolname + git_version, std::string("Added KmerCounts with sdg-workspace add_counts from ")+ws_filename);
    const auto& names(output_ws.get_all_kmer_count_names());
    if (std::find(names.cbegin(), names.cend(), name) != names.cend()) {
        output_ws.kmer_counts.emplace_back(output_ws, counts_filename);
        output_ws.kmer_counts.back().name = name;
        output_ws.dump_to_disk(output + ".bsgws");
    } else {
        throw std::runtime_error(name + " already exists in the workspace");
    }

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
        case WorkspaceFunctions::ADD_DS:
            add_datastore(argc, argv);
            break;
        case WorkspaceFunctions::ADD_COUNTS:
            add_counts(argc, argv);
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
