//
// Created by Bernardo Clavijo (EI) on 03/11/2017.
//

#include <iostream>
#include "PairedReadMapper.hpp"
#include "KMerIDX.h"
#include "FileReader.h"
#include "SMR.h"

#define mem_limit 4
#define GB 1024*1024*1024


/*
 * Pseudo-reader that gets its sequences from the nodes of the graph.
 */

struct GraphNodeReaderParams {
    uint32_t min_length;
    SequenceGraph * sgp;
};

template<typename FileRecord>
class GraphNodeReader {
public:
    explicit GraphNodeReader(GraphNodeReaderParams params, const std::string &filepath) : params(params), numRecords(1) {
        sg=params.sgp;
        min_length=params.min_length;//TODO: use this
    }
    
    bool next_record(FileRecord& rec) {
        if (numRecords<sg->nodes.size()) {
            rec.id = numRecords;
            rec.seq = sg->nodes[numRecords].sequence;
            rec.name = std::to_string(numRecords);
            ++numRecords;
            stats.totalLength += rec.seq.size();
            return true;
        } else return false;
    }
    ReaderStats getSummaryStatistics() {
        stats.totalRecords = numRecords-1;
        return stats;
    }
private:
    SequenceGraph * sg;
    GraphNodeReaderParams params;
    uint32_t numRecords;
    ReaderStats stats;
    uint32_t min_length;
};

/*
 * Mapping of paired end reads.
 * Reads are mapped through a unique k-mer index. R1 and R2 are processed independently.
 * R1 gets odds ids, R2 gets the next even id, so [1,2] and [3,4] are pairs.
 */
void PairedReadMapper::map_reads(std::string filename1, std::string filename2, prmReadType read_type) {
    std::cout<<"Mapping "<<prmReadTypeDesc[read_type]<<" reads from "<<filename1<<" and "<<filename2<<std::endl;
    if (read_type==prmPE){
        std::string pacbio_reads_file(filename1);

        uint64_t maxmem(mem_limit * GB);
        /*
         * Reads a fasta file and generates the set of kmers and which entry number on the fasta they belong to
         * along with their orientation in the form of a int32_t with + for fwd and - for rev, the output is filtered
         * by min < coverage <= max
         */

        const int k =31;
        const int max_coverage=1;
        const std::string output_prefix("smr_files");
        SMR<KmerIDX,
        kmerIDXFactory<FastaRecord>,
        GraphNodeReader<FastaRecord>,
        FastaRecord, GraphNodeReaderParams, KMerIDXFactoryParams> kmerIDX_SMR({1,&sg}, {k}, maxmem, 0, max_coverage,
                                                                          output_prefix);

        std::vector<KmerIDX> unique_kmers;

        // Get the unique_kmers from the file
        unique_kmers = kmerIDX_SMR.read_from_file(output_prefix);

        std::vector<uint64_t > uniqKmer_statistics(kmerIDX_SMR.summaryStatistics());
        std::cout << "Number of " << int(k) << "-kmers seen in assembly " << uniqKmer_statistics[0] << std::endl;
        std::cout << "Number of " << int(k) << "-kmers that appear more than 0 < kmer_coverage <= " << max_coverage
                  << " in assembly " << uniqKmer_statistics[1] << std::endl;
        std::cout << "Number of contigs from the assembly " << uniqKmer_statistics[2] << std::endl;
    }
}