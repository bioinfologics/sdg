#include <iostream>
#include <fstream>
#include <sdglib/datastores/PairedReadsDatastore.hpp>
#include <sdglib/processors/GraphMaker.hpp>
#include <sdglib/workspace/WorkSpace.hpp>
#include <sdglib/batch_counter/BatchKmersCounter.hpp>
#include "cxxopts.hpp"


std::unordered_set<__uint128_t, int128_hash> countKmers(const WorkSpace &ws, int k, int min_coverage, int num_batches) {
    auto kmer_list = BatchKmersCounter::buildKMerCount(k, ws.paired_read_datastores.back(), min_coverage, ".", ".", num_batches);
    std::unordered_set<__uint128_t, int128_hash> kmers;
    for (uint64_t i = 0; i < kmer_list->size; ++i) {
        __uint128_t kmer;
        memcpy(&kmer, &kmer_list->kmers[i], 16);
        kmers.emplace(kmer);
    }
    return kmers;
}

int main(int argc, char * argv[]) {
    std::cout << "Welcome to sdg-dbg" << std::endl << std::endl;
    std::cout << "Git origin: " << GIT_ORIGIN_URL << " -> " << GIT_BRANCH << std::endl;
    std::cout << "Git commit: " << GIT_COMMIT_HASH << std::endl << std::endl;
    std::cout << "Executed command:" << std::endl;
    for (auto i = 0; i < argc; i++) std::cout << argv[i] << " ";
    std::cout << std::endl << std::endl;

    std::string pr_file, lr_file, output_prefix;
    int min_coverage = 5, k = 63, num_batches(0);
    sdglib::OutputLogLevel = sdglib::LogLevels::INFO;
    try {
        cxxopts::Options options("sdg-dbg", "create a DBG graph from a short-read datastore, and populate KCI");

        options.add_options()
                ("help", "Print help")
                ("p,paired_datastore", "input paired read datastore", cxxopts::value<std::string>(pr_file))
                ("b,disk_batches", "number of disk batches to use", cxxopts::value(num_batches))
                //("l,linked_datastore", "input linked read datastore", cxxopts::value<std::string>(lr_filename))
                ("o,output", "output file prefix", cxxopts::value<std::string>(output_prefix));
        options.add_options("Heuristics")
                ("k", "k value for DBG construction (default: 63)", cxxopts::value<int>(k))
                ("c,min_coverage", "minimum coverage for a kmer to include in DBG", cxxopts::value<int>(min_coverage));


        auto result(options.parse(argc, argv));

        if (result.count("help")) {
            std::cout << options.help({"", "Heuristics"}) << std::endl;
            exit(0);
        }

        if (result.count("p") != 1 or result.count("o") != 1) {
            throw cxxopts::OptionException(" please specify input datastore and output prefix");
        }


    } catch (const cxxopts::OptionException &e) {
        std::cout << "Error parsing options: " << e.what() << std::endl << std::endl
                  << "Use option --help to check command line arguments." << std::endl;
        exit(1);
    }

    std::cout << std::endl;

    WorkSpace ws;
    GraphMaker gm2(ws.sdg);
    ws.add_log_entry("Workspace created with sdg-dbg");
    ws.add_log_entry("Origin datastore: " + pr_file);
    ws.paired_read_datastores.emplace_back(ws, pr_file);
    auto kmers2 = countKmers(ws, k, min_coverage, num_batches);
    gm2.new_graph_from_kmerset_trivial128(kmers2, k);
    sdglib::OutputLog() << "DONE! " << ws.sdg.count_active_nodes() << " nodes in graph" << std::endl;
    gm2.tip_clipping(200);
    sdglib::OutputLog() << "DONE! " << ws.sdg.count_active_nodes() << " nodes in graph" << std::endl;
    ws.kci.index_graph();
    ws.kci.start_new_count();
    ws.kci.add_counts_from_datastore(ws.paired_read_datastores.back());
    ws.sdg.write_to_gfa1(output_prefix + "_DBG.gfa");
    ws.dump_to_disk(output_prefix + ".bsgws");
}

