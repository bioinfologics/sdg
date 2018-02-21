//
// Created by Luis Yanes (EI) on 13/11/2017.
//

#include <iostream>
#include <fstream>

#include <sys/stat.h>

#include <sglib/filesystem/check_or_create_directory.h>
#include <sglib/factories/KMerIDXFactory.h>
#include <sglib/readers/FileReader.h>
#include <sglib/SMR.h>
#include "cxxopts.hpp"
#include "sglib/factories/ContigBlockFactory.h"
#include "sglib/SequenceGraph.h"

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
                ("help", "Print help", cxxopts::value<bool>(),"")
                ("r,reference", "reference FASTA", cxxopts::value<std::string>(ref_gfa_filename), "filepath")
                ("l,min_contig_length", "Minimum contig length", cxxopts::value<uint32_t>(min_contig_length)->default_value("1000"),"uint")
                ("m,min_matches", "Minimum kmers to match before calling a block", cxxopts::value<unsigned int>(min_matches)->default_value("1000"),"uint")
                ("o,output", "output directory prefix", cxxopts::value<std::string>(output_prefix), "string")
                ("s,stats_only", "do not dump detailed information", cxxopts::value<bool>(stats_only)->default_value("false"), "uint");

        options.add_options("Performance")
                ("mem_limit", "Memory limit in GB",cxxopts::value<unsigned int>(mem_limit)->default_value("10"), "uint");

        auto result = options.parse(argc, argv);

        if (result.count("help"))
        {
            std::cout << options.help({""}) << std::endl;
            exit(0);
        }

        if (result.count("r")!=1 or result.count("o")!=1) {
            throw cxxopts::OptionException(" please specify input files and output prefix");
        }

    } catch (const cxxopts::OptionException& e)
    {
        std::cout << "Error parsing options: " << e.what() << std::endl << std::endl
                  <<"Use option --help to check command line arguments." << std::endl;
        exit(1);
    }

    if (!sglib::check_or_create_directory(output_prefix)) {
        exit(1);
    }


    std::cout<< "Welcome to gfa-align"<<std::endl<<std::endl;

    /*
     * Load the GFAs
     */
//    SequenceGraph reference_sg, other_sg;
//    reference_sg.load_from_gfa(ref_gfa_filename);
//    other_sg.load_from_gfa(gfa_filename);


    /*
     * Count the unique kmers on both graphs and select the shared unique kmers
     */
    const int k =31;
    const int max_coverage=1;
    SMR<KmerIDX,
    kmerIDXFactory<FastaRecord>,
    FastaReader<FastaRecord>,
    FastaRecord, FastxReaderParams, KMerIDXFactoryParams> ref_kmerIDX_SMR({1}, {k}, {mem_limit*GB, 0,
                                                                          max_coverage,
                                                                          output_prefix});
    std::vector<KmerIDX> vector_unique_kmers(ref_kmerIDX_SMR.read_from_file(ref_gfa_filename));

//    FastaReader<FastaRecord> referenceReader({1}, ref_gfa_filename.substr(0, ref_gfa_filename.length()-4)+".fasta");
//
//
//
//    std::unordered_set<KmerIDX> reference_unique_kmers(std::make_move_iterator(vector_unique_kmers.begin()),
//                                                       std::make_move_iterator(vector_unique_kmers.end()));
//
//    std::cout << "Database generated" << std::endl;
//    GraphNodeReader<FastaRecord> graphReader({1,other_sg}, other_sg.fasta_filename);
//    std::atomic<uint64_t> mapped_count(0),total_count(0);
//    std::cout << "Mapping " << other_sg.filename << " to " << reference_sg.filename << std::endl;
//    std::ofstream ref_to_other_blocks(output_prefix+"ref_to_other.txt");
//#pragma omp parallel shared(graphReader)
//    {
//        ContigBlockFactory<FastaRecord> blockFactory({output_prefix, k, reference_unique_kmers, min_matches});
//        kmerIDXFactory<FastaRecord> kf({k});
//        std::vector<Block> BlockByPos;
//        bool c;
//        FastaRecord read;
//#pragma omp critical(read_record)
//        c = graphReader.next_record(read);
//        while (c) {
//            blockFactory.setFileRecord(read);
//            blockFactory.next_element(BlockByPos);
//            if (!BlockByPos.empty()) {
//                // Sort blocks by number of matches
//                std::sort(BlockByPos.begin(), BlockByPos.end(), Block::byCount());
//                // Report top 5 blocks by number of matches
//#pragma omp critical(report)
//                {
//                    ref_to_other_blocks << read.name << "\t";
//                    std::copy_n(BlockByPos.crbegin(), std::min(5ul, BlockByPos.size()),
//                                std::ostream_iterator<Block>(ref_to_other_blocks, ", "));
//                    ref_to_other_blocks << std::endl;
//                }
//            }
//            auto tc = ++total_count;
//            if (tc % 100000 == 0) std::cout << mapped_count << " / " << tc << std::endl;
//#pragma omp critical(read_record)
//            c = graphReader.next_record(read);
//        }
//    }
}
