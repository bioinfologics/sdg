//
// Created by Bernardo Clavijo (EI) on 03/12/2017.
//

#ifndef BSG_KMERCOMPRESSIONINDEX_HPP
#define BSG_KMERCOMPRESSIONINDEX_HPP

#include <sglib/factories/KMerCountFactory.h>
#include <sglib/readers/SequenceGraphReader.h>
#include <sglib/datastores/PairedReadsDatastore.hpp>
#include "SMR.h"
#include "sglib/graph/SequenceGraph.hpp"

class CStringKMerFactory : protected KMerFactory {
public:
    explicit CStringKMerFactory(uint8_t k) : KMerFactory(k) {};



    ~CStringKMerFactory() {
#pragma omp critical (KMerFactoryDestructor)
        {
            //std::cout << "Bases processed " << bases << "\n";
        }
    }

    // TODO: Adjust for when K is larger than what fits in uint64_t!
    const void create_kmercounts(std::vector<KmerCount> &mers, const char * s) {
        fkmer=0;
        rkmer=0;
        last_unknown=0;
        uint64_t p(0);
        while (s[p]!='\0') {
            //fkmer: grows from the right (LSB)
            //rkmer: grows from the left (MSB)
            fillKBuf(s[p], p, fkmer, rkmer, last_unknown);
            p++;
            if (last_unknown >= K) {
                if (fkmer <= rkmer) {
                    // Is fwd
                    mers.emplace_back(fkmer,1);
                } else {
                    // Is bwd
                    mers.emplace_back(rkmer,1);
                }
            }
        }
    }
    const void create_kmers(std::vector<uint64_t> &mers, const char * s) {
        fkmer=0;
        rkmer=0;
        last_unknown=0;
        uint64_t p(0);
        while (s[p]!='\0') {
            //fkmer: grows from the right (LSB)
            //rkmer: grows from the left (MSB)
            fillKBuf(s[p], p, fkmer, rkmer, last_unknown);
            p++;
            if (last_unknown >= K) {
                if (fkmer <= rkmer) {
                    // Is fwd
                    mers.emplace_back(fkmer);
                } else {
                    // Is bwd
                    mers.emplace_back(rkmer);
                }
            }
        }
    }
    const void create_kmers_direction(std::vector<std::pair<uint64_t,bool>> &mers, const char * s) {
        fkmer=0;
        rkmer=0;
        last_unknown=0;
        uint64_t p(0);
        while (s[p]!='\0') {
            //fkmer: grows from the right (LSB)
            //rkmer: grows from the left (MSB)
            fillKBuf(s[p], p, fkmer, rkmer, last_unknown);
            p++;
            if (last_unknown >= K) {
                if (fkmer <= rkmer) {
                    // Is fwd
                    mers.emplace_back(fkmer,false);
                } else {
                    // Is bwd
                    mers.emplace_back(rkmer,true);
                }
            }
        }
    }
};

class CStringKMerFactory128 : protected KMerFactory128 {
public:
    explicit CStringKMerFactory128(uint8_t k) : KMerFactory128(k) {};



    ~CStringKMerFactory128() {
#pragma omp critical (KMerFactoryDestructor)
        {
            //std::cout << "Bases processed " << bases << "\n";
        }
    }

    const void create_kmers_direction(std::vector<std::pair<__uint128_t,bool>> &mers, const char * s) {
        fkmer=0;
        rkmer=0;
        last_unknown=0;
        uint64_t p(0);
        while (s[p]!='\0') {
            //fkmer: grows from the right (LSB)
            //rkmer: grows from the left (MSB)
            fillKBuf(s[p], p, fkmer, rkmer, last_unknown);
            p++;
            if (last_unknown >= K) {
                if (fkmer <= rkmer) {
                    // Is fwd
                    mers.emplace_back(fkmer,false);
                } else {
                    // Is bwd
                    mers.emplace_back(rkmer,true);
                }
            }
        }
    }
};


class KmerCompressionIndex {
public:
    KmerCompressionIndex(SequenceGraph &_sg, uint64_t _max_mem):sg(_sg){
        max_mem = _max_mem;
    };
    KmerCompressionIndex(SequenceGraph &_sg):sg(_sg){
        max_mem = 0;
    };

    KmerCompressionIndex& operator=(const KmerCompressionIndex &o) {
        if (this != &o) {
            sg = o.sg;
            graph_kmers = o.graph_kmers;
            nodes_depth = o.nodes_depth;
            read_counts = o.read_counts;
            uniq_mode = o.uniq_mode;
        }
        return *this;
    }

    void index_graph();
    void reindex_graph();
    void start_new_count();
    void add_counts_from_file(std::vector<std::string> filename);
    void add_counts_from_datastore(const PairedReadsDatastore & ds);

    void write(std::ofstream & output_file);
    void read(std::ifstream & input_file);
    void save_to_disk(std::string filename);
    void load_from_disk(std::string filename);
    void compute_compression_stats();
    void compute_all_nodes_kci(uint16_t max_graph_freq=10);

    void compute_kci_profiles(std::string filename);

    std::vector<std::vector<uint16_t>> compute_node_coverage_profile(std::string node_sequence, int read_set_index);

    void dump_histogram(std::string filename);

    double compute_compression_for_node(sgNodeID_t node, uint16_t max_graph_freq=10, uint16_t dataset=0);
    SequenceGraph & sg;
    std::vector<KmerCount> graph_kmers;
    std::vector<double> nodes_depth;
    std::vector<std::vector<uint16_t>> read_counts;
    uint64_t max_mem;
    uint16_t uniq_mode=0;

    static const bsgVersion_t min_compat = 0x0001;
};


#endif //BSG_KMERCOMPRESSIONINDEX_HPP
