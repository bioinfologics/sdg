//
// Created by Luis Yanes (EI) on 26/01/2018.
//

#include <iostream>
#include <fstream>
#include <sglib/factories/ContigBlockFactory.h>
#include "cxxopts.hpp"
#include "sglib/SequenceGraph.hpp"
#include "sglib/HaplotypeScorer.hpp"
#include "sglib/Scaffolder.hpp"
#include "sglib/GraphPartitioner.hpp"


int main(int argc, char * argv[]) {
    std::string gfa_filename, bubble_contigs_filename, output_prefix, long_reads;
    std::string dump_mapped, load_mapped;
    std::string cidxreads1, cidxreads2, cidx_dump, cidx_load;
    unsigned int log_level(43);
    uint64_t max_mem_gb(4);
    bool stats_only(false);
    try {
//@formatter:off
        cxxopts::Options options("map-lr", "Long read mapping tool");
        options.add_options()
                ("help", "Print help", cxxopts::value<std::string>(),"")
                ("g,gfa", "input gfa file", cxxopts::value<std::string>(gfa_filename), "filepath")
                ("o,output", "output file prefix", cxxopts::value<std::string>(output_prefix), "path")
                ("log_level", "output log level", cxxopts::value<unsigned int>(log_level), "uint");
        options.add_options("Compression Index Options")
                ("cidxread1", "compression index input reads, left", cxxopts::value<std::string>(cidxreads1), "path")
                ("cidxread2", "compression index input reads, right", cxxopts::value<std::string>(cidxreads2), "path")
                ("dump_cidx", "dump compression index", cxxopts::value<std::string>(cidx_dump),"path")
                ("load_cidx", "load compression index", cxxopts::value<std::string>(cidx_load),"path");
        options.add_options("Long Read Options")
                ("r,long_reads", "input long reads", cxxopts::value<std::string>(long_reads), "filepath")
                ("d,dump_to","dump mapped reads to file",cxxopts::value<std::string>(dump_mapped), "filepath")
                ("l,load_from", "load mapped reads from file", cxxopts::value<std::string>(load_mapped), "filepath")
                ("max_mem", "maximum_memory when mapping (GB, default: 4)", cxxopts::value<uint64_t>(max_mem_gb)->default_value("4"), "GB");
//@formatter:on


        auto result = options.parse(argc, argv);

        if (result.count("help")) {
            std::cout << options.help({""}) << std::endl;
            exit(0);
        }

        if (result.count("g") != 1 or result.count("o") != 1 ) {
            throw cxxopts::OptionException(" please specify input files and output prefix");
        }
        if (result.count("r") != 1 ) {
            throw cxxopts::OptionException(" please specify a long reads file");

        }
        if (result.count("cidxread1") != 1 or result.count("cidxread2") != 1) {
            throw cxxopts::OptionException(" please specify both cidxread1 and cidxread2 files for "
                                                   "the kmer compression indexer");
        }

    } catch (const cxxopts::OptionException &e) {
        std::cout << "Error parsing options: " << e.what() << std::endl << std::endl
                  << "Use option --help to check command line arguments." << std::endl;
        exit(1);
    }
    if (!sglib::check_or_create_directory(output_prefix)) {
        exit(1);
    }

    std::cout << "Welcome to map-lr" << std::endl << std::endl;
    if (gfa_filename.size() <= 4 or gfa_filename.substr(gfa_filename.size() - 4, 4) != ".gfa") {

        throw std::invalid_argument("filename of the gfa input does not end in gfa, it ends in '" +
                                    gfa_filename.substr(gfa_filename.size() - 4, 4) + "'");
    }
    sglib::OutputLogLevel = static_cast<sglib::LogLevels>(log_level);
    max_mem_gb *= GB;
    SequenceGraph sg;
    sg.load_from_gfa(gfa_filename);
    std::vector<PairedReadMapper> mappers;
    mappers.emplace_back(sg);
    mappers.emplace_back(sg);

    SMR<KmerIDX,
            kmerIDXFactory<FastaRecord>,
            GraphNodeReader<FastaRecord>,
            FastaRecord, GraphNodeReaderParams, KMerIDXFactoryParams> kmerIDX_SMR({1, sg}, {31}, {4*GB, 0, 1,
                                                                                                 output_prefix});

    // Get the unique_kmers from the graph into a map
    std::cout << "Indexing graph... " << std::endl;
    std::vector<KmerIDX> uniq_mers(kmerIDX_SMR.process_from_memory());
    std::unordered_set<KmerIDX> uniq_set(uniq_mers.begin(),uniq_mers.end());
    std::unordered_map<uint64_t, graphPosition> kmer_to_graphposition;
    for (auto &kidx :uniq_mers) kmer_to_graphposition[kidx.kmer]={kidx.contigID,kidx.pos};

    /*
     * Map the reads using the smr block factory, transform onto mappings and fill the appropriate vectors on the mapper.
     */

    SMR<Block,
            ContigBlockFactory<FastqRecord>,
            FastqReader<FastqRecord>,
            FastqRecord, FastxReaderParams, FilterSetParams> contigBlock_SMR({0},   // ReaderParams
                                                                            {"out", 31, uniq_set},
                                                                            {4*GB});

    std::vector<Block> read_blocks(contigBlock_SMR.read_from_file(long_reads));


    /*
     * Transform the blocks from vector to mapper
     */

    for (const auto &rb:read_blocks) {
        mappers[0].read_to_node.emplace_back(rb.contigID);
        mappers[0].read_to_node.emplace_back(0);
    }


    // Map long reads or load previous mappings
    if (!load_mapped.empty() and sglib::check_file(load_mapped)) {
        std::cout << "Loading reads from file: " << load_mapped << std::endl;
        mappers[0].load_from_disk(load_mapped);
        mappers[0].print_stats();
    } else {
        std::cout << "Mapping reads... " << std::endl;
        mappers[0].map_reads(long_reads, max_mem_gb);
        if (!dump_mapped.empty()) {
            std::cout << "Saving mappings to file: " << dump_mapped << std::endl;
            mappers[0].save_to_disk(dump_mapped);
        }
    }

    std::cout<<std::endl<<"=== Loading reads compression index ==="<<std::endl;
    //compression index
    KmerCompressionIndex kci(sg,max_mem_gb);
    if (sglib::check_file(cidx_load)){
        kci.load_from_disk(cidx_load);
    } else {
        kci.index_graph();
        kci.start_new_count();
        kci.add_counts_from_file(cidxreads1);
        kci.add_counts_from_file(cidxreads2);
        if (!cidx_dump.empty())
        kci.save_to_disk(cidx_dump);
    }

    if (!kci.read_counts.empty()) {
        kci.compute_compression_stats();
        kci.dump_histogram("kci_histogram.csv");
    }
    // Identify repeaty nodes
    kci.index_graph();
    std::vector<LinkedReadMapper> linkedReadMappers;
    Scaffolder scaf(sg, mappers, linkedReadMappers, kci);
    auto repeatyNodes (scaf.find_repeaty_nodes());
    std::vector<SequenceGraphPath> solvable_paths;
    // For each repeaty node
    for (const auto &central_repeat_node:repeatyNodes) {
        std::cout << "* Evaluating central node " << central_repeat_node << ": " << std::endl;
        // Separate nodes based on reads
        std::vector<sgNodeID_t> fnNodes;
        std::cout << "BW Links" << std::endl;
        for (const auto &l: sg.get_bw_links(central_repeat_node)) {
            std::cout << l << std::endl;
            fnNodes.push_back(l.dest);
        }

        if (std::find(fnNodes.begin(), fnNodes.end(), central_repeat_node) == fnNodes.end())
            fnNodes.push_back(central_repeat_node);

        std::cout << "FW Links" << std::endl;
        for (const auto &l: sg.get_fw_links(central_repeat_node)) {
            std::cout << l << std::endl;
            fnNodes.push_back(l.dest);
        }
        std::vector<sgNodeID_t> fn_repeat(fnNodes.cbegin(), fnNodes.cend());
        std::cout << "Nodes to find partitions for: ";
        std::copy(fn_repeat.cbegin(), fn_repeat.cend(), std::ostream_iterator<sgNodeID_t>(std::cout, " "));
        std::cout << std::endl;

        SequenceSubGraph ssg(sg, fn_repeat);
        GraphPartitioner gp(sg, mappers, kci);
        // Find solutions from the partitions
        auto tp = gp.tags_patterns(ssg);
        auto parts = gp.generate_partitions(ssg, tp, sg.get_bw_links(central_repeat_node).size());
        // Use the solutions to create new paths
        auto parts_score = gp.score_partition_set(ssg, parts, tp);
        auto subgraphs = gp.partitions_as_subgraphs(ssg, parts);
        std::set<sgNodeID_t> tabu_nodes;
        for (auto &psg:subgraphs) {
            std::vector<sgNodeID_t > path_nodes;
            // If there are no three elements (central, bw, fw)
            // Invalid partition
            if (psg.nodes.size() != 3) {
                std::cerr << "Invalid partition size" << std::endl;
                continue;
            } else {
                std::vector<sgNodeID_t > bw_nodes;
                for(const auto &n:sg.get_bw_links(central_repeat_node)) {
                    bw_nodes.push_back(abs(n.dest));
                }
                std::sort(bw_nodes.begin(),bw_nodes.end());
                std::vector<sgNodeID_t > fw_nodes;
                for(const auto &n:sg.get_fw_links(central_repeat_node)) {
                    fw_nodes.push_back(abs(n.dest));
                }
                std::sort(fw_nodes.begin(), fw_nodes.end());
                // check bw nodes
                for (auto &n:psg.nodes) {
                    if (std::find(bw_nodes.begin(), bw_nodes.end(), n) != bw_nodes.end()) path_nodes.push_back(-n);
                    if (std::find(bw_nodes.begin(), bw_nodes.end(), -n) != bw_nodes.end()) path_nodes.push_back(-n);
                }
                if (path_nodes.size() != 1) {
                    std::cerr << "Couldn't find bw node " << std::endl;
                    break;
                }
                path_nodes.push_back(central_repeat_node);
                // check fw nodes
                for (auto &n:psg.nodes) {
                    if (std::find(fw_nodes.begin(), fw_nodes.end(), n) != fw_nodes.end()) path_nodes.push_back(n);
                    if (std::find(fw_nodes.begin(), fw_nodes.end(), -n) != fw_nodes.end()) path_nodes.push_back(n);
                }
                if (path_nodes.size() != 3) {
                    std::cerr << "Couldn't find fw node " << std::endl;
                    break;
                }
            }
            std::cout << "Partitions:" << std::endl;
            unsigned pnumb = 0;
            for (auto &part_node:parts) {
                std::cout << "Partition #" << pnumb++ << ":";
                for (auto n:part_node) std::cout << " " << (n ? 1 : 0);
                std::cout << std::endl;
            }
            std::cout << "Partition solutions as nodes:" << std::endl;
            pnumb = 0;
            std::cout << "Partition #" << pnumb++ << ":";
            for (auto n:psg.nodes) std::cout << " " << n;
            std::cout << std::endl;

            std::cout << "Nodes to link ";
            std::copy(path_nodes.begin(),path_nodes.end(), std::ostream_iterator<sgNodeID_t>(std::cout, ", "));
            std::cout << std::endl;
            std::cout << "Links: " << std::endl;
            for (const auto &n:path_nodes){
                auto fw_links(sg.get_fw_links(n));
                auto bw_links(sg.get_bw_links(n));
                std::copy(fw_links.begin(),fw_links.end(), std::ostream_iterator<Link>(std::cout, "\n"));
                std::copy(bw_links.begin(),bw_links.end(), std::ostream_iterator<Link>(std::cout, "\n"));
            }
            // Construct a valid path through the nodes.
            // Join valid path using the SequenceGraph helper.
            bool tabued(false);
            for (const auto &n: path_nodes)
                if (tabu_nodes.find(n) != tabu_nodes.end()) {tabued=true;break;}
            if (tabued) break;

            // Add solvable path to the list of paths
            solvable_paths.emplace_back(sg,path_nodes);
//            sg.join_path(SequenceGraphPath(sg, path_nodes));
            tabu_nodes.insert(path_nodes.begin(), path_nodes.end());
        }
//        std::cout << "Scores: " << parts_score.first << " " << parts_score.second << std::endl;
    }
    std::cout << repeatyNodes.size() << " repeat nodes generated " << solvable_paths.size() << " solvable paths" << std::endl;

    sg.write_to_gfa("salida.gfa");
}
