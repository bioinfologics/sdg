//
// Created by Luis Yanes (EI) on 26/01/2018.
//

#include <iostream>
#include <fstream>
#include <sglib/factories/ContigBlockFactory.h>
#include <sglib/filesystem/helpers.h>
#include <sglib/readers/FileReader.h>
#include <sglib/readers/SequenceGraphReader.h>
#include <sglib/SMR.h>
#include <sglib/mappers/LongReadMapper.hpp>
#include "cxxopts.hpp"

#include <sglib/utilities/omp_safe.hpp>
#include <vector>

void map_using_lrMapper(uint8_t k, uint8_t w, SequenceGraph &sg, std::string &output_prefix, std::string &long_reads) {
    if (k > 28) {
        throw (std::invalid_argument(std::string("Unsupported K=") + std::to_string(k) + std::string(" value, K must be lower then 28")));
    }
    (w==0)?(uint8_t)(k*.66f):w;
    LongReadsDatastore datastore(long_reads);
    LongReadMapper rm(sg, datastore, k);
    rm.update_graph_index();
    rm.map_reads();
}

int main(int argc, char * argv[]) {
    sglib::OutputLog(false) << "Welcome to LongRead Mapper"<<std::endl<<std::endl;
    sglib::OutputLog(false) << "Git origin: " << GIT_ORIGIN_URL << " -> "  << GIT_BRANCH << std::endl;
    sglib::OutputLog(false) << "Git commit: " << GIT_COMMIT_HASH << std::endl<<std::endl;
    sglib::OutputLog() << "Executed command:"<<std::endl;
    for (auto i=0;i<argc;i++) sglib::OutputLog(false) <<argv[i]<<" ";
    sglib::OutputLog(false) <<std::endl<<std::endl;

    std::string gfa_filename, bubble_contigs_filename, output_prefix, lr_datastore;
    std::string dump_mapped, load_mapped;
    unsigned int log_level(0);
    uint64_t max_mem_gb(4);
    bool stats_only(false);
    uint8_t K(15), W(0);
//@formatter:off
    cxxopts::Options options("map-lr", "LongRead Mapper");
    options.add_options()
            ("help", "Print help", cxxopts::value<std::string>(),"")
            ("k,mer_size", "K-mer size for indexing/mapping", cxxopts::value<uint8_t>(K)->default_value("15"), "0-28")
            ("w,window_size", "Minimiser window size defaults to 2/3 of K value", cxxopts::value<uint8_t>(W), "0-19")
            ("g,gfa", "input gfa file", cxxopts::value<std::string>(gfa_filename), "filepath")
            ("o,output", "output file prefix", cxxopts::value<std::string>(output_prefix), "path")
            ("l,log_level", "output log level", cxxopts::value<unsigned int>(log_level)->default_value("0"), "uint")
            ("d,datastore", "long read datastore", cxxopts::value<std::string>(lr_datastore), "filepath")
            ("max_mem", "maximum_memory when mapping (GB, default: 4)", cxxopts::value<uint64_t>(max_mem_gb)->default_value("4"), "GB");
//@formatter:on
    try {
        auto result = options.parse(argc, argv);

        if (result.count("help")) {
            std::cout << options.help({""}) << std::endl;
            exit(0);
        }
        if (K > 28) {
            throw cxxopts::OptionException(std::string("Unsupported K=") + std::to_string(K) + std::string(" value, K must be lower then 28"));
        }
        (W==0)?(uint8_t)(K*.66f):W;

        if (result.count("g") != 1 or result.count("o") != 1 ) {
            throw cxxopts::OptionException(" please specify input files and output prefix");
        }
    } catch (const cxxopts::OptionException &e) {
        std::cout << "Error parsing options: " << e.what() << std::endl;
        std::cout << options.help({""}) << std::endl;
        exit(1);
    }
    if (!sglib::check_or_create_directory(output_prefix)) {
        exit(1);
    }

    if (gfa_filename.size() <= 4 or gfa_filename.substr(gfa_filename.size() - 4, 4) != ".gfa") {

        throw std::invalid_argument("filename of the gfa input does not end in gfa, it ends in '" +
                                    gfa_filename.substr(gfa_filename.size() - 4, 4) + "'");
    }
    sglib::OutputLogLevel = static_cast<sglib::LogLevels>(log_level);
    max_mem_gb *= GB;
    SequenceGraph sg;
    sg.load_from_gfa(gfa_filename);
    LongReadsDatastore datastore(lr_datastore);
    LongReadMapper rm(sg, datastore, K);

    rm.update_graph_index();
    rm.map_reads();

    return 0;
}
