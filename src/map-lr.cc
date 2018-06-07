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
#include <sglib/mappers/minimap2/minimap.h>

void
printMatch(const mm_idx_t *mi, std::ofstream &matchOutput, uint32_t readID, const std::string &read_name,
           int read_len, const mm_reg1_t *regs0, int j) {
    matchOutput << readID
                << "\t" << read_name
                << "\t" << read_len
                << "\t" << j + 1
                << "\t" << regs0[j].rid
                << "\t" << mi->seq[regs0[j].rid].name
                << "\t" << mi->seq[regs0[j].rid].len
                << "\t" << regs0[j].rs
                << "\t" << regs0[j].re
                << "\t" << regs0[j].qs
                << "\t" << regs0[j].qe
                << "\t" << "+-"[regs0[j].rev]
                << "\t" << regs0[j].mlen
                << "\t" << regs0[j].blen
                << "\t" << regs0[j].mapq
                << "\n";
}

void map_reads(SequenceGraph &sg, mm_mapopt_t *opt, mm_idx_t *mi, LongReadsDatastore &datastore, std::unordered_set<uint32_t> readIDs) {
    {
        std::ofstream matchHeader("header.paf");
        matchHeader << "read_id\t"
                    << "read_name\t"
                    << "read_len\t"
                    << "read_nmap\t"
                    << "node_id\t"
                    << "node_name\t"
                    << "node_len\t"
                    << "node_start\t"
                    << "node_end\t"
                    << "read_start\t"
                    << "read_end\t"
                    << "strand\t"
                    << "nmatches\t"
                    << "blocklen\t"
                    << "mapq\n";
    }
#pragma omp parallel
{
    mm_tbuf_t *buf = mm_tbuf_init();
    std::ofstream matchOutput(std::string("thread_")+std::to_string(omp_get_thread_num())+std::string(".paf"));
#pragma omp for
        for (uint32_t readID = 1; readID < datastore.size(); ++readID) {
            std::string read_seq(datastore.get_read_sequence(readID));
            std::string read_name(std::to_string(readID));
            int read_len(static_cast<int>(read_seq.length()));
            int n_regs0;
            mm_reg1_t *regs0 = mm_map(mi, read_len, read_seq.data(), &n_regs0, buf, opt, read_name.data());

            if (n_regs0<=1) {
                for (int j = 0; j < n_regs0; j++) {
                    printMatch(mi, matchOutput, readID, read_name, read_len, regs0, j);
                }
            }
            if (n_regs0 > 1) {
                for (int j = 0; j < n_regs0 - 1; ++j) {
                    auto fromNode = sg.oldnames_to_ids.at(mi->seq[regs0[j].rid].name);
                    auto toNode = sg.oldnames_to_ids.at(mi->seq[regs0[j+1].rid].name);

                    if (std::abs(fromNode) != std::abs(toNode)) {
                        if (!sg.link_exists(fromNode, toNode)) {
                            continue;
                        }
                    }
                    printMatch(mi, matchOutput, readID, read_name, read_len, regs0, j);
                    printMatch(mi, matchOutput, readID, read_name, read_len, regs0, j+1);
                }
            }
            // Check that the links exists in the graph.
            // Store the path in a path list

            // Report if a read has broken paths (and store the missing links).


            for (int i = 0; i<n_regs0;i++) free(regs0[i].p);
            free(regs0);
        }
    mm_tbuf_destroy(buf);
}

}

void map_using_minimap(uint8_t k, SequenceGraph &sg, std::string &lr_datastore) {
    std::vector<const char *> names(sg.nodes.size());
    std::vector<const char *> seqs(sg.nodes.size());
    for (std::vector<std::string>::size_type i = 1; i < seqs.size(); i++) {
        names[i] = sg.oldnames[i].data();
        seqs[i] = sg.nodes[i].sequence.data();
    }

    mm_mapopt_t opt;
    mm_mapopt_init(&opt);
    mm_idx_t *graph_index = mm_idx_str(15, k, 0, 14, static_cast<int>(seqs.size()-1), &seqs[1], &names[1]);
    mm_mapopt_update(&opt, graph_index);
    opt.flag |= MM_F_CIGAR;

    LongReadsDatastore datastore(lr_datastore);
    map_reads(sg, &opt, graph_index, datastore, {});
    mm_idx_destroy(graph_index);
}

void map_using_lrMapper(uint8_t k, uint8_t w, SequenceGraph &sg, std::string &output_prefix, std::string &long_reads) {
    if (k > 28) {
        throw (std::invalid_argument(std::string("Unsupported K=") + std::to_string(k) + std::string(" value, K must be lower then 28")));
    }
    (w==0)?(uint8_t)(k*.66f):w;
    LongReadsDatastore datastore(long_reads);
    LongReadMapper rm(k, w, sg, datastore);
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
    LongReadMapper rm(K, W, sg, datastore);


    if (0) {
        map_using_minimap(K, sg, lr_datastore);
    }
    if (1) {
        rm.update_graph_index();
        rm.map_reads();
    }

    return 0;
}
