//
// Created by Bernardo Clavijo (EI) on 03/12/2017.
//

#ifndef BSG_KMERCOMPRESSIONINDEX_HPP
#define BSG_KMERCOMPRESSIONINDEX_HPP

#include <sglib/factories/KMerCountFactory.h>
#include <sglib/readers/SequenceGraphReader.h>
#include "SMR.h"
#include "SequenceGraph.h"

class KmerCompressionIndex {
public:
    KmerCompressionIndex(SequenceGraph &_sg, uint64_t _max_mem):sg(_sg){
        max_mem = _max_mem;
    };
    KmerCompressionIndex(SequenceGraph &_sg):sg(_sg){
        max_mem = 0;
    };
    void index_graph();
    void reindex_graph();
    void start_new_count();
    void add_counts_from_file(std::vector<std::string> filename);

    void write(std::ofstream & output_file);
    void read(std::ifstream & input_file);
    void save_to_disk(std::string filename);
    void load_from_disk(std::string filename);
    void compute_compression_stats();
    void compute_all_nodes_kci(uint16_t max_graph_freq=10);

    void dump_histogram(std::string filename);

    double compute_compression_for_node(sgNodeID_t node, uint16_t max_graph_freq=10, uint16_t dataset=0);
    SequenceGraph & sg;
    std::vector<KmerCount> graph_kmers;
    std::vector<double> nodes_depth;
    std::vector<std::vector<uint16_t>> read_counts;
    uint64_t max_mem;
    uint16_t uniq_mode=0;

};


#endif //BSG_KMERCOMPRESSIONINDEX_HPP
