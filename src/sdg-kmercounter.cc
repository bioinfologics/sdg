#include <iostream>
#include <fstream>
#include <sdglib/graph/SequenceDistanceGraph.hpp>
#include <sdglib/workspace/WorkSpace.hpp>
#include "cxxopts.hpp"

const std::string program_name("sdg-kmercounts");
enum class KmerCountsFunctions{
    NONE, //This is needed because the default of a map is 0
    MAKE,
    ADD,
    STATS,
};
struct KmerCountsFunctionMap : public std::map<std::string, KmerCountsFunctions>
{
    KmerCountsFunctionMap()
    {
        operator[]("make") =  KmerCountsFunctions::MAKE;
        operator[]("add") =  KmerCountsFunctions::ADD;
        operator[]("stats") = KmerCountsFunctions::STATS;
    };
};

void add_kmer_count(int argc, char **argv) {
    int k(31);
    std::vector<std::string> fastq_files;
    std::string output;
    std::string ds_filename;
    std::string counts_filename;
    std::string name;
    try {
        cxxopts::Options options(program_name + " add", "SDG add a KmerCount to a KmerCounter file");

        options.add_options()
                ("n,name", "KmerCounter name", cxxopts::value(name))
                ("c,counts_file", "input counts file (used for index)", cxxopts::value(counts_filename))
                ("f,fastq", "input reads (multi)", cxxopts::value(fastq_files))
                ("d,datastore", "input datastore", cxxopts::value(ds_filename))
                ("k", "kmer length", cxxopts::value(k)->default_value("31"))
                ("o,output", "output counts prefix", cxxopts::value(output))
                ("help", "Print help");
        auto newargc=argc-1;
        auto newargv=&argv[1];
        auto result=options.parse(newargc,newargv);
        if (result.count("help")) {
            std::cout << options.help({""}) << std::endl;
            exit(0);
        }

        if (output.empty()) {
            throw cxxopts::OptionException(" please specify an output prefix for the counts file");
        }
        if (name.empty()) {
            throw cxxopts::OptionException(" please specify a name for this count");
        }

        if (counts_filename.empty()) {
            throw cxxopts::OptionException(" please specify an input counts file");
        }

        if ( (result.count("fastq")<1 and ds_filename.empty()) or (!result.count("fastq")<1 and !ds_filename.empty())) {
            throw cxxopts::OptionException(" please specify an input files (datastore or fastqs)");
        }

    } catch (const cxxopts::OptionException &e) {
        std::cout << "Error parsing options: " << e.what() << std::endl << std::endl
                  << "Use option --help to check command line arguments." << std::endl;
        exit(1);
    }

    WorkSpace ws;
    KmerCounter kc(ws, counts_filename);
    if (!fastq_files.empty()) {
        kc.add_count(name, fastq_files);
    }
    if (!ds_filename.empty()) {
        if (ds_filename.substr(ds_filename.find(".") + 1) == "prseq") {
            kc.add_count(name, PairedReadsDatastore(ws, ds_filename));
        }
        if (ds_filename.substr(ds_filename.find(".") + 1) == "lrseq") {
            kc.add_count(name, LinkedReadsDatastore(ws, ds_filename));
        }
        if (ds_filename.substr(ds_filename.find(".") + 1) == "loseq") {
            kc.add_count(name, LongReadsDatastore(ws, ds_filename));
        }
    }
    std::ofstream output_file(output+".sdgkc");
    kc.write_counts(output_file);

}

void make_kmer_counts(int argc, char **argv) {
    int k(27);
    std::vector<std::string> fastq_files;
    std::string output;
    std::string gfa_filename;
    std::string ws_filename;
    std::string ds_filename;
    std::string counts_filename;
    std::string name;
    try {
        cxxopts::Options options(program_name + " make", "SDG make KmerCounter");

        options.add_options()
                ("w,workspace", "input workspace (used for index)", cxxopts::value(ws_filename))
                ("g,gfa", "input gfa file (used for index)", cxxopts::value(gfa_filename))
                ("n,name", "KmerCounter name", cxxopts::value(name))
                ("c,counts_file", "input counts file (used for index)", cxxopts::value(counts_filename))
                ("f,fastq", "input reads (multi)", cxxopts::value(fastq_files))
                ("d,datastore", "input datastore", cxxopts::value(ds_filename))
                ("k", "kmer length", cxxopts::value(k)->default_value("27"))
                ("o,output", "output counts prefix", cxxopts::value(output))
                ("help", "Print help");
        auto newargc=argc-1;
        auto newargv=&argv[1];
        auto result=options.parse(newargc,newargv);
        if (result.count("help")) {
            std::cout << options.help({""}) << std::endl;
            exit(0);
        }

        if (output.empty()) {
            throw cxxopts::OptionException(" please specify an output prefix for the counts file");
        }
        if (name.empty()) {
            throw cxxopts::OptionException(" please specify a name for this count");
        }

        if ( (result.count("fastq")<1 and ds_filename.empty()) or (!result.count("fastq")<1 and !ds_filename.empty())) {
            throw cxxopts::OptionException(" please specify an input files (datastore or fastqs)");
        }

        int base_count=0;
        gfa_filename.empty() ? base_count : base_count++;
        ws_filename.empty() ? base_count : base_count++;
        counts_filename.empty() ? base_count : base_count++;
        if ( base_count== 0 or base_count > 1 ) {
            throw cxxopts::OptionException(" please specify a gfa, workspace or counts file, this file is required for the index");
        }
    } catch (const cxxopts::OptionException &e) {
        std::cout << "Error parsing options: " << e.what() << std::endl << std::endl
                  << "Use option --help to check command line arguments." << std::endl;
        exit(1);
    }

    //===== LOAD GRAPH =====
    if (!ws_filename.empty() or !gfa_filename.empty()) {
        WorkSpace ws;
        if (!ws_filename.empty()) {
            ws.load_from_disk(ws_filename);
        } else if (!gfa_filename.empty()) {
            SequenceDistanceGraph sg(ws);
            sg.load_from_gfa(gfa_filename);
        }
        KmerCounter kc(ws, name, k);
        if (!fastq_files.empty()) {
            kc.add_count(name, fastq_files);
        }
        if (!ds_filename.empty()) {
            if (ds_filename.substr(ds_filename.find(".") + 1) == "prseq") {
                kc.add_count(name, PairedReadsDatastore(ws, ds_filename));
            }
            if (ds_filename.substr(ds_filename.find(".") + 1) == "lrseq") {
                kc.add_count(name, LinkedReadsDatastore(ws, ds_filename));
            }
            if (ds_filename.substr(ds_filename.find(".") + 1) == "loseq") {
                kc.add_count(name, LongReadsDatastore(ws, ds_filename));
            }
        }
        std::ofstream output_file(output+".sdgkc");
        kc.write_counts(output_file);
    } else if (!counts_filename.empty()) {
        WorkSpace ws;
        KmerCounter kc(ws, counts_filename);
        if (!fastq_files.empty()) {
            kc.add_count(name, fastq_files);
        }
        if (!ds_filename.empty()) {
            if (ds_filename.substr(ds_filename.find(".") + 1) == "prseq") {
                kc.add_count(name, PairedReadsDatastore(ws, ds_filename));
            }
            if (ds_filename.substr(ds_filename.find(".") + 1) == "lrseq") {
                kc.add_count(name, LinkedReadsDatastore(ws, ds_filename));
            }
            if (ds_filename.substr(ds_filename.find(".") + 1) == "loseq") {
                kc.add_count(name, LongReadsDatastore(ws, ds_filename));
            }
        }
        std::ofstream output_file(output+".sdgkc");
        kc.write_counts(output_file);
    }
}

void stats_kmer_counts(int argc, char **argv) {
    std::string filename;
    try {

        cxxopts::Options options(program_name +" stats", "SDG KmerCounter stats");

        options.add_options()
                ("help", "Print help")
                ("n,name", "KmerCounter name", cxxopts::value<std::string>(filename));

        auto newargc=argc-1;
        auto newargv=&argv[1];
        auto result=options.parse(newargc,newargv);
        if (result.count("help")) {
            std::cout << options.help({""}) << std::endl;
            exit(0);
        }

        if (result.count("name")==0) {
            throw cxxopts::OptionException(" please specify KmerCounter file");
        }


    } catch (const cxxopts::OptionException &e) {
        std::cout << "Error parsing options: " << e.what() << std::endl << std::endl
                  << "Use option --help to check command line arguments." << std::endl;
        exit(1);
    }
    WorkSpace ws;

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
        case KmerCountsFunctions::ADD:
            add_kmer_count(argc, argv);
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

