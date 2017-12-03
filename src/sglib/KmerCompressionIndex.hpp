//
// Created by Bernardo Clavijo (EI) on 03/12/2017.
//

#ifndef BSG_KMERCOMPRESSIONINDEX_HPP
#define BSG_KMERCOMPRESSIONINDEX_HPP

#include <sglib/factories/KMerCountFactory.h>
#include <sglib/readers/SequenceGraphReader.h>
#include "SMR.h"
#include "SequenceGraph.hpp"

class KmerCompressionIndex {
public:
    KmerCompressionIndex(SequenceGraph &_sg, uint64_t max_mem);
    void index_graph();
    void start_new_count();
    void add_counts_from_file(std::string filename);

    void save_to_disk(std::string filename);
    void load_from_disk(std::string filename);
    void compute_compression_stats();

    void compute_compression_for_node(sgNodeID_t node);
    SequenceGraph & sg;
    std::vector<KmerCount> graph_kmers;
    std::vector<std::vector<uint16_t>> read_counts;
    uint64_t max_mem;
    uint16_t uniq_mode=0;

};


#endif //BSG_KMERCOMPRESSIONINDEX_HPP
