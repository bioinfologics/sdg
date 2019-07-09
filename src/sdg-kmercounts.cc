#include <iostream>
#include <fstream>
#include <sdglib/processors/KmerCompressionIndex.hpp>
#include <sdglib/graph/SequenceDistanceGraph.hpp>
#include <sdglib/workspace/WorkSpace.hpp>
#include "cxxopts.hpp"

const std::string program_name("sdg-kmerscounts");
enum class KmerCountsFunctions{
    NONE, //This is needed because the default of a map is 0
    MAKE,
    STATS,
};
struct KmerCountsFunctionMap : public std::map<std::string, KmerCountsFunctions>
{
    KmerCountsFunctionMap()
    {
        operator[]("make") =  KmerCountsFunctions::MAKE;
        operator[]("stats") = KmerCountsFunctions::STATS;
    };
};

void make_kmer_counts(int argc, char **argv) {
    std::vector<std::string> fastq_files;
    std::string output;
    std::string gfa_filename;
    std::string ws_filename;
    std::string name;
    try {
        cxxopts::Options options(program_name + " make", "SDG make KmerCounts");

        options.add_options()
                ("help", "Print help")
                ("n,name", "KmerCounts name", cxxopts::value(name))
                ("w,workspace", "input workspace", cxxopts::value<std::string>(ws_filename))
                ("g,gfa", "input gfa file", cxxopts::value<std::string>(gfa_filename))
                ("f,fastq", "input reads (multi)", cxxopts::value<std::vector<std::string>>(fastq_files))
                ("o,output", "output file", cxxopts::value<std::string>(output));
        auto newargc=argc-1;
        auto newargv=&argv[1];
        auto result=options.parse(newargc,newargv);
        if (result.count("help")) {
            std::cout << options.help({""}) << std::endl;
            exit(0);
        }

        if (result.count("fastq")<1 or output=="") {
            throw cxxopts::OptionException(" please specify input files and output prefix");
        }

        if (gfa_filename.empty() and ws_filename.empty()) {
            throw cxxopts::OptionException(" please specify a gfa, or workspace");
        }

    } catch (const cxxopts::OptionException &e) {
        std::cout << "Error parsing options: " << e.what() << std::endl << std::endl
                  << "Use option --help to check command line arguments." << std::endl;
        exit(1);
    }

    //===== LOAD GRAPH =====
    WorkSpace ws;
    if (!ws_filename.empty()) {
        ws.load_from_disk(ws_filename);
        ws.add_kmer_counts_datastore()
    } else if (!gfa_filename.empty()) {

        SequenceDistanceGraph sg(ws);
    }
}

void stats_kmer_counts(int argc, char **argv) {
    std::string filename;
    try {

        cxxopts::Options options(program_name +" stats", "SDG KmerCounts stats");

        options.add_options()
                ("help", "Print help")
                ("n,name", "KmerCounts name", cxxopts::value<std::string>(filename));

        auto newargc=argc-1;
        auto newargv=&argv[1];
        auto result=options.parse(newargc,newargv);
        if (result.count("help")) {
            std::cout << options.help({""}) << std::endl;
            exit(0);
        }

        if (result.count("name")==0) {
            throw cxxopts::OptionException(" please specify KmerCounts file");
        }


    } catch (const cxxopts::OptionException &e) {
        std::cout << "Error parsing options: " << e.what() << std::endl << std::endl
                  << "Use option --help to check command line arguments." << std::endl;
        exit(1);
    }
    WorkSpace ws;
    SequenceDistanceGraph sg(ws);

}

int main(int argc, char * argv[]) {

    std::cout << program_name <<std::endl<<std::endl;
    std::cout << "Git origin: " << GIT_ORIGIN_URL << " -> "  << GIT_BRANCH << std::endl;
    std::cout << "Git commit: " << GIT_COMMIT_HASH << std::endl<<std::endl;
    std::cout << "Executed command:"<<std::endl;

    KmerCountsFunctionMap functionMap;

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
        case KmerCountsFunctions::MAKE:
            make_kmer_counts(argc, argv);
            break;
        case KmerCountsFunctions::STATS:
            stats_kmer_counts(argc, argv);
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

