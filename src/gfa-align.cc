//
// Created by Luis Yanes (EI) on 13/11/2017.
//

#include <iostream>
#include <fstream>
#include <algorithm>
#include <sglib/KMerIDX.h>
#include <sglib/FileReader.h>
#include <sglib/SMR.h>
#include "cxxopts.hpp"
#include "sglib/ContigBlocks.h"
#include "sglib/SequenceGraph.hpp"

int main(int argc, char * argv[]) {
    std::string gfa_filename,ref_gfa_filename,output_prefix;
    std::uint32_t min_contig_length(0);
    unsigned int mem_limit(10);
    unsigned int min_matches(1);
    bool stats_only=0;
    try
    {
        cxxopts::Options options("gfa-align", "GFA Alignment tool");

        options.add_options()
                ("help", "Print help")
                ("g,gfa", "input GFA", cxxopts::value<std::string>(gfa_filename), "filepath")
                ("r,reference", "reference GFA", cxxopts::value<std::string>(ref_gfa_filename), "filepath")
                ("l,min_contig_length", "Minimum contig length", cxxopts::value<uint32_t>(min_contig_length)->default_value("1000"),"uint")
                ("m,min_matches", "Minimum kmers to match before calling a block", cxxopts::value<unsigned int>(min_matches)->default_value("1000"),"uint")
                ("o,output", "output file prefix", cxxopts::value<std::string>(output_prefix), "string")
                ("s,stats_only", "do not dump detailed information", cxxopts::value<bool>(stats_only)->default_value("0"), "uint");

        options.add_options("Performance")
                ("mem_limit", "Memory limit in GB",cxxopts::value<unsigned int>(mem_limit)->default_value("10"), "uint");
        options.parse(argc, argv);

        if (options.count("help"))
        {
            std::cout << options.help({""}) << std::endl;
            exit(0);
        }

        if (options.count("g")!=1 or options.count("r")!=1 or options.count("o")!=1) {
            throw cxxopts::OptionException(" please specify input files and output prefix");
        }

    } catch (const cxxopts::OptionException& e)
    {
        std::cout << "Error parsing options: " << e.what() << std::endl << std::endl
                  <<"Use option --help to check command line arguments." << std::endl;
        exit(1);
    }


    std::cout<< "Welcome to gfa-align"<<std::endl<<std::endl;

    /*
     * Load the GFAs
     */
    SequenceGraph reference_sg, other_sg;
    reference_sg.load_from_gfa(ref_gfa_filename);
    other_sg.load_from_gfa(gfa_filename);


    /*
     * Count the unique kmers on both graphs and select the shared unique kmers
     */
    const int k =31;
    const int max_coverage=1;
    const std::string smr_output_prefix("smr_files");
    SMR<KmerIDX,
    kmerIDXFactory<FastaRecord>,
    GraphNodeReader<FastaRecord>,
    FastaRecord, GraphNodeReaderParams, KMerIDXFactoryParams> ref_kmerIDX_SMR({1,reference_sg}, {k}, mem_limit*GB, 0, max_coverage,
                                                                          smr_output_prefix+"_ref");

    std::vector<KmerIDX> reference_unique_kmers;

    // Get the unique_kmers from the file
    reference_unique_kmers = ref_kmerIDX_SMR.read_from_file(output_prefix);

    GraphNodeReader<FastaRecord> graphReader({1,other_sg}, other_sg.fasta_filename);
    std::atomic<uint64_t> mapped_count(0),total_count(0);
    ContigBlockFactory<FastaRecord> blockFactory({output_prefix, k, reference_unique_kmers, min_matches});
    std::vector<std::vector<Block>> ctg_ctgMatches;
#pragma omp parallel shared(graphReader,blockFactory)
    {
        std::vector<KmerIDX> readkmers;
        kmerIDXFactory<FastaRecord> kf({k});
        std::vector<Block> BlockByPos;
        bool c;
        FastaRecord read;
#pragma omp critical(read_record)
        c = graphReader.next_record(read);
        while (c) {
            blockFactory.setFileRecord(read);
            blockFactory.next_element(BlockByPos);
            if (!BlockByPos.empty()) {
                std::cout << read.name << "\t";
                std::sort(BlockByPos.begin(), BlockByPos.end(), Block::byCount());
                std::copy_n(BlockByPos.cbegin(), std::min(5ul,BlockByPos.size()), std::ostream_iterator<Block>(std::cout, ", "));
                std::cout << std::endl;
            }
#pragma omp critical(add_blocks)
            ctg_ctgMatches.push_back(BlockByPos);
            auto tc = ++total_count;
            if (tc % 100000 == 0) std::cout << mapped_count << " / " << tc << std::endl;
#pragma omp critical(read_record)
            c = graphReader.next_record(read);
        }
    }

//    for (auto it = ctg_ctgMatches.cbegin(); it != ctg_ctgMatches.cend(); ++it){
//        if (it->empty()) continue;
//        std::cout << std::distance(ctg_ctgMatches.cbegin(), it) << "\t";
//        std::cout << std::endl;
//    }
}
