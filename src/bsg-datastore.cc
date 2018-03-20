#include <iostream>
#include <fstream>
#include <sglib/datastores/LinkedReadsDatastore.hpp>
#include "cxxopts.hpp"


int main(int argc, char * argv[]) {
    std::cout << "bsg-datastore"<<std::endl<<std::endl;
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
        std::string read1, read2, read_type, output;
        uint16_t readsize=150;
        try {
            cxxopts::Options options("bsg-datastore make", "BSG make datastore");

            options.add_options()
                    ("help", "Print help")
                    ("1,read1", "input reads, left", cxxopts::value<std::string>(read1))
                    ("2,read2", "input reads, right", cxxopts::value<std::string>(read2))
                    ("t,read_type", "One of: paired,10x", cxxopts::value<std::string>(read_type))
                    ("s,max_read_size", "max read size for short reads, truncates if longer (default 150)", cxxopts::value<uint16_t>(readsize))
                    ("o,output", "output file", cxxopts::value<std::string>(output));
            auto newargc=argc-1;
            auto newargv=&argv[1];
            auto result=options.parse(newargc,newargv);
            if (result.count("help")) {
                std::cout << options.help({""}) << std::endl;
                exit(0);
            }

            if (read1 == "" or read2 == "" or read_type == "" or output == "") {
                throw cxxopts::OptionException(" please specify input files, type and output prefix");
            }


        } catch (const cxxopts::OptionException &e) {
            std::cout << "Error parsing options: " << e.what() << std::endl << std::endl
                      << "Use option --help to check command line arguments." << std::endl;
            exit(1);
        }

        //===== DATASTORE CREATION =====
        if (read_type == "10x" or read_type == "10xseq") {
            LinkedReadsDatastore ds(read1, read2, output+".lrseq",(read_type == "10xseq" ? LinkedReadsFormat::seq
                                                                                 : LinkedReadsFormat::UCDavis),readsize);
            //ds.dump_index_to_disk(output+".lrIdx");
        } else {
            std::cout << "read_type '" << read_type << "' is not supported (yet?)" << std::endl;
        }

    }
    else if (0==strcmp(argv[1],"view")) {
        std::vector<std::string> filenames;
        try {

            cxxopts::Options options("bsg-datastore view", "BSG view datastore");

            options.add_options()
                    ("help", "Print help")
                    ("d,datastore", "datastore name (multi)", cxxopts::value<std::vector<std::string>>(filenames));

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
        std::cout<<"Datastore updates not implemented yet."<<std::endl;
    }
    else {
        std::cout<<"Please specify one of: make, update, view"<<std::endl;
    }
}

