#include <iostream>
#include <fstream>
#include <sglib/KmerCompressionIndex.hpp>
#include <sglib/WorkSpace.hpp>
#include "cxxopts.hpp"


int main(int argc, char * argv[]) {
    std::cout << "bsg-workspace"<<std::endl<<std::endl;
    std::cout << "Git origin: " << GIT_ORIGIN_URL << " -> "  << GIT_BRANCH << std::endl;
    std::cout << "Git commit: " << GIT_COMMIT_HASH << std::endl<<std::endl;
    std::cout << "Executed command:"<<std::endl;
    for (auto i=0;i<argc;i++) std::cout<<argv[i]<<" ";
    std::cout<<std::endl<<std::endl;

    if (argc <2){
        std::cout<<"Please specify one of: make, log, status, dump"<<std::endl;
        exit(1);
    }

    if (0==strcmp(argv[1],"make")) {
        std::vector<std::string> lr_datastores;
        std::string output="";
        std::string gfa_filename="",kci_filename="";
        try {
            cxxopts::Options options("bsg-kmerspectra make", "BSG make workspace");

            options.add_options()
                    ("help", "Print help")
                    ("g,gfa", "input gfa file", cxxopts::value<std::string>(gfa_filename))
                    ("l,linked_reads", "linked reads datastore", cxxopts::value<std::vector<std::string>>(lr_datastores))
                    ("k,kmerspectra", "KCI kmer spectra", cxxopts::value<std::string>(kci_filename))
                    ("o,output", "output file", cxxopts::value<std::string>(output));
            auto newargc=argc-1;
            auto newargv=&argv[1];
            auto result=options.parse(newargc,newargv);
            if (result.count("help")) {
                std::cout << options.help({""}) << std::endl;
                exit(0);
            }

            if (output=="" or gfa_filename=="") {
                throw cxxopts::OptionException(" please specify gfa and output prefix");
            }


        } catch (const cxxopts::OptionException &e) {
            std::cout << "Error parsing options: " << e.what() << std::endl << std::endl
                      << "Use option --help to check command line arguments." << std::endl;
            exit(1);
        }

        //===== LOAD GRAPH =====
        WorkSpace w;
        w.add_log_entry("Created with bsg-makeworkspace");
        w.sg.load_from_gfa(gfa_filename);
        w.add_log_entry("GFA imported from "+gfa_filename+" ("+std::to_string(w.sg.nodes.size()-1)+" nodes)");
        if (kci_filename!=""){
            w.kci.load_from_disk(kci_filename);
            w.add_log_entry("KCI kmer spectra imported from "+kci_filename);
        }
        for (auto lrds:lr_datastores){
            //create and load the datastore, and the mapper!
            w.linked_read_datastores.emplace_back(lrds);
            w.linked_read_mappers.emplace_back(w.sg,w.linked_read_datastores.back());
            w.add_log_entry("LinkedReadDatastore imported from "+lrds+" ("+std::to_string(w.linked_read_datastores.back().size())+" reads)");
        }


        w.dump_to_disk(output);

    }
    else if (0==strcmp(argv[1],"log")) {
        std::string filename;
        try {

            cxxopts::Options options("bsg-workspace log", "BSG workspace log");

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
                throw cxxopts::OptionException(" please specify kmer spectra file");
            }


        } catch (const cxxopts::OptionException &e) {
            std::cout << "Error parsing options: " << e.what() << std::endl << std::endl
                      << "Use option --help to check command line arguments." << std::endl;
            exit(1);
        }
        WorkSpace w;
        w.load_from_disk(filename,true);
        w.print_log();
    }
    else if (0==strcmp(argv[1],"dump")) {
        std::string filename,gfafilename,nodeinfofilename;
        try {

            cxxopts::Options options("bsg-workspace log", "BSG workspace log");

            options.add_options()
                    ("help", "Print help")
                    ("w,workspace", "workspace filename", cxxopts::value<std::string>(filename))
                    ("g,gfa", "gfa output prefix", cxxopts::value<std::string>(gfafilename))
                    ("n,node_info", "node info prefix",cxxopts::value<std::string>(nodeinfofilename));

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
        WorkSpace w;
        w.load_from_disk(filename);
        if (!w.sg.is_sane()) {
            sglib::OutputLog()<<"ERROR: sg.is_sane() = false"<<std::endl;
            //return 1;
        }
        if (not gfafilename.empty()){
            w.sg.write_to_gfa(gfafilename+".gfa");
        }
        if (not nodeinfofilename.empty()){
           std::ofstream nif(nodeinfofilename+".csv");
           nif<<"ID, lenght, kci"<<std::endl;
           for (auto n=1;n<w.sg.nodes.size();++n){
               if (w.sg.nodes[n].status==sgNodeStatus_t::sgNodeDeleted) continue;
               nif<<n<<", "<<w.sg.nodes[n].sequence.size()<<", "<<w.kci.compute_compression_for_node(n,1)<<std::endl;
           }
        }
    }
    else if (0==strcmp(argv[1],"node-kci-dump")){
        std::string filename;
        sgNodeID_t cnode;
        std::string prefix;

        try {
            cxxopts::Options options("bsg-workspace kci-dump", "BSG workspace kci-dump");

            options.add_options()
                    ("help", "Print help")
                    ("w,workspace", "workspace filename", cxxopts::value<std::string>(filename))
                    ("n,node", "node to profile", cxxopts::value<sgNodeID_t>(cnode))
                    ("p,prefix", "Prefix for the output file", cxxopts::value<std::string>(prefix));

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

        std::cout << "Here" << filename << std::endl;

        WorkSpace w;
        w.load_from_disk(filename);
        
        std::ofstream reads_ofl("reads_profile"+std::to_string(cnode)+".cvg", std::ios_base::app);
        std::ofstream unique_ofl("unique_profile"+std::to_string(cnode)+".cvg", std::ios_base::app);
        std::ofstream assm_ofl("assmcn_profile"+std::to_string(cnode)+".cvg", std::ios_base::app);

        std::cout << "Profile para el nodo:" << cnode << std::endl;
        std::string sequence = w.sg.nodes[cnode].sequence;

        auto coverages = w.kci.compute_node_coverage_profile(sequence);

        std::cout << ">Reads_"<< prefix << "_" <<cnode << std::endl;
        reads_ofl << ">Reads_"<< prefix << "_"<< cnode << "|";
        for (auto c: coverages[0]){
            std::cout << c << " ";
            reads_ofl << c << " ";
        }
        std::cout << std::endl;
        reads_ofl << std::endl;

        std::cout << ">Uniqueness_"<< prefix << "_"<<cnode << std::endl;
        unique_ofl << ">Uniqueness_"<< prefix << "_" << cnode << "|";
        for (auto c: coverages[1]){
            std::cout << c << " ";
            unique_ofl << c << " ";
        }
        std::cout << std::endl;
        unique_ofl << std::endl;

        std::cout << ">Graph_"<< prefix << "_" << cnode << std::endl;
        assm_ofl << ">Graph_"<< prefix << "_"<<cnode << "|";
        for (auto c: coverages[2]){
            std::cout << c << " ";
            assm_ofl << c << " ";
        }
        std::cout << std::endl;
        assm_ofl << std::endl;

        reads_ofl.close();
        unique_ofl.close();
        assm_ofl.close();
    }
    else if (0==strcmp(argv[1], "kci-stack")) {
        // Take multiple kci files and stack them in the reads kmer vector for later use

    }
    else {
        std::cout<<"Please specify one of: make, log, status, dump"<<std::endl;
    }
}

