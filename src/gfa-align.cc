//
// Created by Luis Yanes (EI) on 13/11/2017.
//

#include <iostream>
#include <fstream>
#include <sglib/KMerIDX.h>
#include <sglib/FileReader.h>
#include <sglib/SMR.h>
#include "cxxopts.hpp"
#include "sglib/ContigBlocks.h"
#include "sglib/SequenceGraph.hpp"

int main(int argc, char * argv[]) {
    std::string gfa_filename,ref_gfa_filename,output_prefix;
    std::uint32_t min_contig_length(0);
    bool stats_only=0;

    try
    {
        cxxopts::Options options("gfa-align", "GFA Alignment tool");

        options.add_options()
                ("help", "Print help")
                ("g,gfa", "input GFA", cxxopts::value<std::string>(gfa_filename), "filepath")
                ("r,reference", "reference GFA", cxxopts::value<std::string>(ref_gfa_filename), "filepath")
                ("m,min_contig_length", "Minimum contig length", cxxopts::value<uint32_t>(min_contig_length)->default_value("1000"),"uint")
                ("o,output", "output file prefix", cxxopts::value<std::string>(output_prefix), "string")
                ("s,stats_only", "do not dump detailed information", cxxopts::value<bool>(stats_only)->default_value("0"));

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
    reference_sg.load_from_gfa(ref_gfa_filename+".gfa");
    other_sg.load_from_gfa(gfa_filename+".gfa");


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

    /*
 * Store the lengths of the sequences longer than min_contig_length
 * and keep a set of the shorter ones to filter the unique_kmers map with
 */
    FastaReader<FastaRecord> reader({0}, gfa_filename+".fasta");
    std::unordered_map<int32_t, std::pair<uint64_t, uint32_t>> sequences_fasta;
    std::unordered_set<int> invalid_sequences;
    FastaRecord rec;
    while (reader.next_record(rec)) {
        if (rec.seq.size() > min_contig_length) {
            sequences_fasta[rec.id] = std::make_pair(rec.seq.size(),0);
        } else {
            invalid_sequences.insert(rec.id);
        }
    }

    /*
     * Delete all the kmers in the map that come from "invalid" unitigs
     */
    std::remove_if(total_uniq_set.begin(), total_uniq_set.end(),
                   [&invalid_sequences] (const KmerIDX &kidx)
                   {
                       // Remove if I found the contigID on the invalid_sequences
                       return (invalid_sequences.find(kidx.contigID) != invalid_sequences.cend());
                   }
    );
    std::cout << "Number of " << int(k) << "-mers in contigs larger than " << min_contig_length << " " << total_uniq_set.size() << std::endl;
    // Cleanup invalid_sequences
    invalid_sequences = std::unordered_set<int>();

    for (const auto &kdx: total_uniq_set) {
        sequences_fasta[std::abs(kdx.contigID)].second++;
    }


    FastaReader<FastaRecord> fastaReader({min_contig_length}, gfa_filename+".fasta");
    std::atomic<uint64_t> mapped_count(0),total_count(0);
    ContigBlockFactory<FastaRecord> blockFactory({output_prefix, k, total_uniq_set, sequences_fasta});
    std::vector<std::vector<Block>> total_validBlocks;
#pragma omp parallel shared(fastaReader,blockFactory)
    {
        std::vector<KmerIDX> readkmers;
        kmerIDXFactory<FastaRecord> kf({k});
        std::vector<Block> validBlocks;
        bool c;
        FastaRecord read;
#pragma omp critical(read_record)
        c = fastaReader.next_record(read);
        while (c) {
            blockFactory.setFileRecord(read);
            blockFactory.next_element(validBlocks);
#pragma omp critical(add_blocks)
            total_validBlocks.push_back(validBlocks);
            auto tc=++total_count;
            if (tc % 100000 == 0) std::cout << mapped_count << " / " << tc << std::endl;
#pragma omp critical(read_record)
            c = fastaReader.next_record(read);
        }
    }



    return 0;
}
