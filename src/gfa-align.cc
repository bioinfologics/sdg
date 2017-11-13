//
// Created by Luis Yanes (EI) on 13/11/2017.
//

#include <iostream>
#include <fstream>
#include <sglib/KMerIDX.h>
#include <sglib/FileReader.h>
#include <sglib/SMR.h>
#include "cxxopts.hpp"
#include "sglib/SequenceGraph.hpp"

int main(int argc, char * argv[]) {
    std::string gfa_filename,ref_gfa_filename,output_prefix;
    bool stats_only=0;

    try
    {
        cxxopts::Options options("gfa-align", "GFA Alignment tool");

        options.add_options()
                ("help", "Print help")
                ("g,gfa", "input GFA", cxxopts::value<std::string>(gfa_filename))
                ("r,reference", "reference GFA", cxxopts::value<std::string>(ref_gfa_filename))
                ("o,output", "output file prefix", cxxopts::value<std::string>(output_prefix))
                ("s,stats_only", "do not dump detailed information (default=0)", cxxopts::value<bool>(stats_only));

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
    const int min_matches=1;
    const std::string smr_output_prefix("smr_files");
    SMR<KmerIDX,
    kmerIDXFactory<FastaRecord>,
    GraphNodeReader<FastaRecord>,
    FastaRecord, GraphNodeReaderParams, KMerIDXFactoryParams> ref_kmerIDX_SMR({1,&reference_sg}, {k}, 4*GB, 0, max_coverage,
                                                                          smr_output_prefix+"_ref");

    std::vector<KmerIDX> reference_unique_kmers;

    // Get the unique_kmers from the file
    reference_unique_kmers = ref_kmerIDX_SMR.read_from_file(output_prefix);

    SMR<KmerIDX,
            kmerIDXFactory<FastaRecord>,
            GraphNodeReader<FastaRecord>,
            FastaRecord, GraphNodeReaderParams, KMerIDXFactoryParams> asm_kmerIDX_SMR({1,&reference_sg}, {k}, 4*GB, 0, max_coverage,
                                                                                  smr_output_prefix+"_asm");

    std::vector<KmerIDX> asm_unique_kmers;

    // Get the unique_kmers from the file
    asm_unique_kmers = asm_kmerIDX_SMR.read_from_file(output_prefix);

    std::vector<KmerIDX> total_uniq_set;
    std::set_intersection(reference_unique_kmers.cbegin(),reference_unique_kmers.cend(), asm_unique_kmers.cbegin(), asm_unique_kmers.cend(), std::back_inserter(total_uniq_set));

    std::cout << reference_unique_kmers.size() << " " << asm_unique_kmers.size() << " " << total_uniq_set.size() << std::endl;


    return 0;
}
