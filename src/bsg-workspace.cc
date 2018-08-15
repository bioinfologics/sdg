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
        std::vector<std::string> lr_datastores,pr_datastores,Lr_datastores;
        std::string output="";
        std::string gfa_filename="",kci_filename="";
        try {
            cxxopts::Options options("bsg-kmerspectra make", "BSG make workspace");

            options.add_options()
                    ("help", "Print help")
                    ("g,gfa", "input gfa file", cxxopts::value<std::string>(gfa_filename))
                    ("p,paired_reads", "paired reads datastore", cxxopts::value<std::vector<std::string>>(pr_datastores))
                    ("l,linked_reads", "linked reads datastore", cxxopts::value<std::vector<std::string>>(lr_datastores))
                    ("L,long_reads", "long reads datastore", cxxopts::value<std::vector<std::string>>(Lr_datastores))
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
        for (auto prds:pr_datastores){
            //create and load the datastore, and the mapper!
            w.paired_read_datastores.emplace_back(prds);
            w.paired_read_mappers.emplace_back(w.sg,w.paired_read_datastores.back());
            w.add_log_entry("PairedReadDatastore imported from "+prds+" ("+std::to_string(w.paired_read_datastores.back().size())+" reads)");
        }
        for (auto lrds:lr_datastores){
            //create and load the datastore, and the mapper!
            w.linked_read_datastores.emplace_back(lrds);
            w.linked_read_mappers.emplace_back(w.sg,w.linked_read_datastores.back());
            w.add_log_entry("LinkedReadDatastore imported from "+lrds+" ("+std::to_string(w.linked_read_datastores.back().size())+" reads)");
        }
        for (auto Lrds:Lr_datastores){
            //create and load the datastore, and the mapper!
            w.long_read_datastores.emplace_back(Lrds);
            w.long_read_mappers.emplace_back(w.sg,w.long_read_datastores.back());
            w.add_log_entry("LongReadDatastore imported from "+Lrds+" ("+std::to_string(w.long_read_datastores.back().size())+" reads)");
        }


        w.dump_to_disk(output+".bsgws");

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
        std::string filename,gfafilename,nodeinfofilename,seqfilename;
        float minKCI=.5, maxKCI=1.25;
        size_t min_size=2000,max_size=100000000;
        try {

            cxxopts::Options options("bsg-workspace log", "BSG workspace log");

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
        if (not seqfilename.empty()) {
            std::ofstream sof(seqfilename + ".fasta");
            for (auto n:w.select_from_all_nodes(min_size,max_size,0,UINT32_MAX, minKCI, maxKCI)){
                sof<<">seq"<<n<<std::endl<<w.sg.nodes[n].sequence<<std::endl;
            }
        }

    }
    else if (0==strcmp(argv[1],"copy_mapped")) {
        std::vector<int> lr_datastores,pr_datastores,Lr_datastores;
        std::string output="";
        std::string base_filename="",merge_filename="";
        bool force(false);
        try {
            cxxopts::Options options("bsg-kmerspectra make", "BSG make workspace");

            options.add_options()
                    ("help", "Print help")
                    ("w,workspace_base", "base workspace", cxxopts::value<std::string>(base_filename))
                    ("m,workspace_merge", "merge workspace", cxxopts::value<std::string>(merge_filename))
                    ("p,paired_reads", "paired reads datastore", cxxopts::value<std::vector<int>>(pr_datastores))
                    ("l,linked_reads", "linked reads datastore", cxxopts::value<std::vector<int>>(lr_datastores))
                    ("L,long_reads", "long reads datastore", cxxopts::value<std::vector<int>>(Lr_datastores))
                    ("f,force", "force copy of datastores", cxxopts::value<bool>(force))
                    ("o,output", "output file", cxxopts::value<std::string>(output));
            auto newargc=argc-1;
            auto newargv=&argv[1];
            auto result=options.parse(newargc,newargv);
            if (result.count("help")) {
                std::cout << options.help({""}) << std::endl;
                exit(0);
            }

            if (output=="" or base_filename=="" or merge_filename == "") {
                throw cxxopts::OptionException(" please specify a base workspace, a merge workspace and output prefix");
            }

        } catch (const cxxopts::OptionException &e) {
            std::cout << "Error parsing options: " << e.what() << std::endl << std::endl
                      << "Use option --help to check command line arguments." << std::endl;
            exit(1);
        }

        WorkSpace base,merge,out;
        base.load_from_disk(base_filename);
        merge.load_from_disk(merge_filename);
        out.kci = base.kci;
        out.sg = base.sg;
        std::set<int> lr_filter(lr_datastores.begin(), lr_datastores.end());
        std::set<int> pr_filter(pr_datastores.begin(), pr_datastores.end());
        std::set<int> Lr_filter(Lr_datastores.begin(), Lr_datastores.end());

        std::unordered_set<std::string> base_datastores;
        for (int i = 0; i < base.paired_read_datastores.size(); ++i) {
            base_datastores.insert(base.paired_read_datastores[i].filename);
            if (lr_filter.empty() or lr_filter.count(i) == 0) {
                sglib::OutputLog()<< "Adding " << merge.paired_read_datastores[i].filename << std::endl;
                out.paired_read_datastores.push_back(base.paired_read_datastores[i]);
                out.paired_read_mappers.push_back(base.paired_read_mappers[i]);
            }
        }

        for (int i = 0; i < base.linked_read_datastores.size(); ++i) {
            base_datastores.insert(base.linked_read_datastores[i].filename);
            if (lr_filter.empty() or lr_filter.count(i) == 0) {
                sglib::OutputLog()<< "Adding " << base.linked_read_datastores[i].filename << std::endl;
                out.linked_read_datastores.push_back(base.linked_read_datastores[i]);
                out.linked_read_mappers.push_back(base.linked_read_mappers[i]);
            }
        }

        for (int i = 0; i < base.long_read_datastores.size(); ++i) {
            base_datastores.insert(base.long_read_datastores[i].filename);
            if (lr_filter.empty() or lr_filter.count(i) == 0) {
                sglib::OutputLog()<< "Adding " << base.long_read_datastores[i].filename << std::endl;
                out.long_read_datastores.push_back(base.long_read_datastores[i]);
                out.long_read_mappers.push_back(base.long_read_mappers[i]);
            }
        }

        if (base.sg == merge.sg or force) {
            for (int i = 0; i < merge.paired_read_datastores.size(); ++i) {
                if (base_datastores.find(merge.paired_read_datastores[i].filename) == base_datastores.end()) {
                    sglib::OutputLog()<< "Adding " << merge.paired_read_datastores[i].filename << std::endl;
                    out.paired_read_datastores.push_back(merge.paired_read_datastores[i]);
                    out.paired_read_mappers.push_back(merge.paired_read_mappers[i]);
                }
            }
            for (int i = 0; i < merge.linked_read_datastores.size(); ++i) {
                if (base_datastores.find(merge.linked_read_datastores[i].filename) == base_datastores.end()) {
                    sglib::OutputLog()<< "Adding " << merge.linked_read_datastores[i].filename << std::endl;
                    out.linked_read_datastores.push_back(merge.linked_read_datastores[i]);
                    out.linked_read_mappers.push_back(merge.linked_read_mappers[i]);
                }
            }
            for (int i = 0; i < merge.long_read_datastores.size(); ++i) {
                if (base_datastores.find(merge.long_read_datastores[i].filename) == base_datastores.end()) {
                    sglib::OutputLog()<< "Adding " << merge.long_read_datastores[i].filename << std::endl;
                    out.long_read_datastores.push_back(merge.long_read_datastores[i]);
                    out.long_read_mappers.push_back(merge.long_read_mappers[i]);
                }
            }
        } else {
            sglib::OutputLog() << "The graphs in the base and merge datastores are different, not merging mappers"
                               << std::endl;
            for (int i = 0; i < merge.paired_read_datastores.size(); ++i) {
                if (base_datastores.find(merge.paired_read_datastores[i].filename) == base_datastores.end()) {
                    sglib::OutputLog()<< "Adding " << merge.paired_read_datastores[i].filename << " without mappings" << std::endl;
                    out.paired_read_datastores.push_back(merge.paired_read_datastores[i]);
                    out.paired_read_mappers.emplace_back(merge.sg,merge.paired_read_datastores[i]);
                }
            }
            for (int i = 0; i < merge.linked_read_datastores.size(); ++i) {
                if (base_datastores.find(merge.linked_read_datastores[i].filename) == base_datastores.end()) {
                    sglib::OutputLog()<< "Adding " << merge.linked_read_datastores[i].filename << " without mappings" << std::endl;
                    out.linked_read_datastores.push_back(merge.linked_read_datastores[i]);
                    out.linked_read_mappers.emplace_back(merge.sg,merge.linked_read_datastores[i]);
                }
            }
            for (int i = 0; i < merge.long_read_datastores.size(); ++i) {
                if (base_datastores.find(merge.long_read_datastores[i].filename) == base_datastores.end()) {
                    sglib::OutputLog()<< "Adding " << merge.long_read_datastores[i].filename << " without mappings" << std::endl;
                    out.long_read_datastores.push_back(merge.long_read_datastores[i]);
                    out.long_read_mappers.emplace_back(merge.sg,merge.long_read_datastores[i]);
                }
            }
        }
        out.dump_to_disk(output+".bsgws");
    }
    else {
        std::cout<<"Please specify one of: make, log, status, dump, copy_mapped"<<std::endl;
    }
}

