#include <iostream>
#include <fstream>
#include <sglib/KmerCompressionIndex.hpp>
#include "cxxopts.hpp"


int main(int argc, char * argv[]) {
    std::cout << "bsg-kmerspectra"<<std::endl<<std::endl;
    std::cout << "Git origin: " << GIT_ORIGIN_URL << " -> "  << GIT_BRANCH << std::endl;
    std::cout << "Git commit: " << GIT_COMMIT_HASH << std::endl<<std::endl;
    std::cout << "Executed command:"<<std::endl;
    for (auto i=0;i<argc;i++) std::cout<<argv[i]<<" ";
    std::cout<<std::endl<<std::endl;

    if (argc <2){
        std::cout<<"Please specify one of: make, update, view"<<std::endl;
        exit(1);
    }

    if (0==strcmp(argv[1],"make")) {
        std::vector<std::string> fastq_files;
        std::string output;
        std::string gfa_filename;
        try {
            cxxopts::Options options("bsg-kmerspectra make", "BSG make kmers pectra");

            options.add_options()
                    ("help", "Print help")
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

            if (result.count("fastq")<1 or output=="" or gfa_filename=="") {
                throw cxxopts::OptionException(" please specify gfa, input files, output prefix");
            }


        } catch (const cxxopts::OptionException &e) {
            std::cout << "Error parsing options: " << e.what() << std::endl << std::endl
                      << "Use option --help to check command line arguments." << std::endl;
            exit(1);
        }

        //===== LOAD GRAPH =====
        SequenceGraph sg;
        KmerCompressionIndex kci(sg);
        sg.load_from_gfa(gfa_filename);
        kci.index_graph();
        kci.start_new_count();
        kci.add_counts_from_file(fastq_files);
        kci.save_to_disk(output);



    }
    else if (0==strcmp(argv[1],"view")) {
        std::vector<std::string> filenames;
        try {

            cxxopts::Options options("bsg-kmerspectra view", "BSG view kmer spectra");

            options.add_options()
                    ("help", "Print help")
                    ("s,kmerspectra", "kmerspectra name (multi)", cxxopts::value<std::vector<std::string>>(filenames));

            auto newargc=argc-1;
            auto newargv=&argv[1];
            auto result=options.parse(newargc,newargv);
            if (result.count("help")) {
                std::cout << options.help({""}) << std::endl;
                exit(0);
            }

            if (result.count("datastore")==0) {
                throw cxxopts::OptionException(" please specify datastore file (s)");
            }


        } catch (const cxxopts::OptionException &e) {
            std::cout << "Error parsing options: " << e.what() << std::endl << std::endl
                      << "Use option --help to check command line arguments." << std::endl;
            exit(1);
        }



    }
    else if(0==strcmp(argv[1],"update")){
        std::cout<<"Kmerspectra updates not implemented yet."<<std::endl;
    }
    else {
        std::cout<<"Please specify one of: make, update, view"<<std::endl;
    }
}

