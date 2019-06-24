#include <iostream>
#include <fstream>
#include <sdglib/processors/KmerCompressionIndex.hpp>
#include "cxxopts.hpp"


int main(int argc, char * argv[]) {
    std::cout << "sdg-kmerspectra"<<std::endl<<std::endl;
    std::cout << "Git origin: " << GIT_ORIGIN_URL << " -> "  << GIT_BRANCH << std::endl;
    std::cout << "Git commit: " << GIT_COMMIT_HASH << std::endl<<std::endl;
    std::cout << "Executed command:"<<std::endl;
    for (auto i=0;i<argc;i++) std::cout<<argv[i]<<" ";
    std::cout<<std::endl<<std::endl;

    if (argc <2){
        std::cout<<"Please specify one of: make, stats, apply, transform"<<std::endl;
        exit(1);
    }

    if (0==strcmp(argv[1],"make")) {
        std::vector<std::string> fastq_files;
        std::string output;
        std::string gfa_filename;
        try {
            cxxopts::Options options("sdg-kmerspectra make", "BSG make kmers pectra");

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
        SequenceDistanceGraph sg;
        KmerCompressionIndex kci(sg);
        sg.load_from_gfa(gfa_filename);
        kci.index_graph();
        kci.start_new_count();
        kci.add_counts_from_file(fastq_files);
        kci.compute_compression_stats();
        kci.save_to_disk(output);



    }
    else if (0==strcmp(argv[1],"stats")) {
        std::string filename;
        try {

            cxxopts::Options options("sdg-kmerspectra stats", "BSG kmer spectra stats");

            options.add_options()
                    ("help", "Print help")
                    ("s,kmerspectra", "kmerspectra name", cxxopts::value<std::string>(filename));

            auto newargc=argc-1;
            auto newargv=&argv[1];
            auto result=options.parse(newargc,newargv);
            if (result.count("help")) {
                std::cout << options.help({""}) << std::endl;
                exit(0);
            }

            if (result.count("kmerspectra")==0) {
                throw cxxopts::OptionException(" please specify kmer spectra file");
            }


        } catch (const cxxopts::OptionException &e) {
            std::cout << "Error parsing options: " << e.what() << std::endl << std::endl
                      << "Use option --help to check command line arguments." << std::endl;
            exit(1);
        }
        SequenceDistanceGraph sg;
        KmerCompressionIndex kci(sg);
        kci.load_from_disk(filename);
        sdglib::OutputLog()<<kci.graph_kmers.size()<<" kmers in spectra loaded from disk"<<std::endl;
        kci.compute_compression_stats();


    }
    else if(0==strcmp(argv[1],"apply")){
        std::string kci_filename,gfa_filename,output;
        int maxfreq=10;
        bool reindex=false;
        try {

            cxxopts::Options options("sdg-kmerspectra apply", "BSG kmer spectra apply");

            options.add_options()
                    ("help", "Print help")
                    ("g,gfa", "input gfa file", cxxopts::value<std::string>(gfa_filename))
                    ("s,kmerspectra", "kmerspectra name", cxxopts::value<std::string>(kci_filename))
                    ("f,max_freq", "maximum graph frequency (default:10)", cxxopts::value<int>(maxfreq))
                    ("r,reindex", "re-index the graph counts (default: false)", cxxopts::value<bool>(reindex))
                    ("o,output", "output gfa file with depth", cxxopts::value<std::string>(output));

            auto newargc=argc-1;
            auto newargv=&argv[1];
            auto result=options.parse(newargc,newargv);
            if (result.count("help")) {
                std::cout << options.help({""}) << std::endl;
                exit(0);
            }

            if (result.count("kmerspectra")==0 or result.count("gfa")==0 or result.count("output")==0) {
                throw cxxopts::OptionException(" please specify gfa, kmer spectra and output files");
            }


        } catch (const cxxopts::OptionException &e) {
            std::cout << "Error parsing options: " << e.what() << std::endl << std::endl
                      << "Use option --help to check command line arguments." << std::endl;
            exit(1);
        }
        SequenceDistanceGraph sdg;
        KmerCompressionIndex kci(sdg);
        sdg.load_from_gfa(gfa_filename);
        kci.load_from_disk(kci_filename);
        sdglib::OutputLog()<<kci.graph_kmers.size()<<" kmers in spectra loaded from disk"<<std::endl;
        if (reindex) kci.reindex_graph();
        kci.compute_compression_stats();
        kci.compute_all_nodes_kci(maxfreq);
        sdg.write_to_gfa1(output, {}, kci.nodes_depth);
    }
    else if (0==strcmp(argv[1],"transform")) {
        std::string kci_filename,gfa_filename,output;

        try {

            cxxopts::Options options("sdg-kmerspectra transform", "BSG kmer spectra transform");

            options.add_options()
                    ("help", "Print help")
                    ("s,kmerspectra", "kmerspectra name", cxxopts::value<std::string>(kci_filename))
                    ("o,output", "output gfa file with depth", cxxopts::value<std::string>(output));

            auto newargc=argc-1;
            auto newargv=&argv[1];
            auto result=options.parse(newargc,newargv);
            if (result.count("help")) {
                std::cout << options.help({""}) << std::endl;
                exit(0);
            }

            if (result.count("kmerspectra")==0 or result.count("output")==0) {
                throw cxxopts::OptionException(" please specify gfa, kmer spectra and output files");
            }


        } catch (const cxxopts::OptionException &e) {
            std::cout << "Error parsing options: " << e.what() << std::endl << std::endl
                      << "Use option --help to check command line arguments." << std::endl;
            exit(1);
        }


        std::ifstream kci_file(kci_filename);
        SequenceDistanceGraph sg;
        KmerCompressionIndex kci(sg);
        uint64_t kcount;
        kci_file.read(( char *) &kcount,sizeof(kcount));
        kci.graph_kmers.resize(kcount);
        kci_file.read(( char *) kci.graph_kmers.data(),sizeof(KmerCount)*kcount);
        //read-to-node
        uint64_t ccount;
        kci_file.read(( char *) &ccount,sizeof(ccount));
        for (auto i=0;i<ccount;++i) {
            kci.read_counts.emplace_back();
            kci.read_counts.back().resize(kcount);
            kci_file.read(( char *) kci.read_counts.back().data(), sizeof(uint16_t) * kcount);
        }
        kci.compute_compression_stats();
        kci.save_to_disk(output);
    }
    else {
        std::cout<<"Please specify one of: make, stats, view, transform"<<std::endl;
    }
}

