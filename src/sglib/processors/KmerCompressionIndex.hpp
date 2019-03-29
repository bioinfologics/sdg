//
// Created by Bernardo Clavijo (EI) on 03/12/2017.
//

#ifndef BSG_KMERCOMPRESSIONINDEX_HPP
#define BSG_KMERCOMPRESSIONINDEX_HPP

#include <sglib/factories/KMerCountFactory.hpp>
#include <sglib/readers/SequenceGraphReader.hpp>
#include <sglib/datastores/PairedReadsDatastore.hpp>
#include "sglib/SMR.hpp"
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
            fillKBuf(s[p], fkmer, rkmer, last_unknown);
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
            fillKBuf(s[p], fkmer, rkmer, last_unknown);
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
            fillKBuf(s[p], fkmer, rkmer, last_unknown);
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
            fillKBuf(s[p], fkmer, rkmer, last_unknown);
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

    /**
     * Adds a new empty kmer count collection at the back of read_counts
     */
    void start_new_count();

    /**
     * Accumulates the kmer count from the provided fastq file to the last available read_counts collection
     * @param path to fastq file
     */
    void add_counts_from_file(std::vector<std::string> filename);

    /**
     * Accumulates the kmer count from the provided data-store to the last available read_counts collection
     * @param ds PairedReadsDatastore ds
     */
    void add_counts_from_datastore(const PairedReadsDatastore & ds);

    void write(std::ofstream & output_file);
    void read(std::ifstream & input_file);
    void save_to_disk(std::string filename);
    void load_from_disk(std::string filename);

    /**
     *  Reports the mean, median and mode for the for the read_counts[0] collection
     */
    void compute_compression_stats();

    /**
     * Computes the kci for each node using the compute_compression_for_node function, the kci is calculated
     * using the read_counts[0] collection. The result is stored in nodes_depth.
     * Only the kmers with graph frequency < max_graph_freq are used in the calculation
     * @param max_graph_freq Maximum kci represented in the graph
     */
    void compute_all_nodes_kci(uint16_t max_graph_freq=10);

    /**
     * Computes the absent kmers profile for all the non deleted nodes of the graph in each available read_count collection
     * For each node executes the function compute_node_coverage_profile in every available read_count collection and
     * returns the number of kmers of each node in each collection considered absent (freq<3)
     * The output is a file with the sequence name and the counts in the same order as they were loaded in the workspace
     * @param filename filename for the output file
     */
    void compute_kci_profiles(std::string filename);

    /**
     * Computes the kmer coverage of the sequence using the read_count[read_set_index] collection and the graph_kmers collection
     * returns a vector of triplets [[read_count, unique_kmers_count, graph_count],...] one for each kmer of the sequence
     * @param node_sequence sequence to compute the coverage
     * @param read_set_index index of the desired read_count collection
     * @return vector of triplets [read_count, unique_kmers_count, graph_count]
     */
    std::vector<std::vector<uint16_t>> compute_node_coverage_profile(std::string node_sequence, int read_set_index);

    /**
     * Dumps the kmer coverage histogram for the graph kmers using the count in read_count[0] collection
     * @param filename
     */
    void dump_histogram(std::string filename);

    /**
     * Computes the kci index for the selected node using the kmers from the node with kmer count in the
     * graph < max_graph_freq, the kmer count in read_count[dataset] collection and the frequency for the main heterozygous peak
     * stored in uniq_mode.
     * This function will ignore the ovelapping kmers
     * if the length og the node without overlap is < k return nan ((double)0/0)
     * @param node node if for the selecte dnode
     * @param max_graph_freq maximum frequency of a kmer in the graph to be considered
     * @param dataset index of the read_count collection to use
     * @return kci for the selected node
     */
    double compute_compression_for_node(sgNodeID_t node, uint16_t max_graph_freq=10, uint16_t dataset=0); // Only for dataset 0

    void print_status();

    SequenceGraph & sg;

    /**
     * Collection of graph kmers
     */
    std::vector<KmerCount> graph_kmers;

    /**
     * computed kci for all nodes in the graph, this vector is filled using compute_all_nodes_kci
     */
    std::vector<double> nodes_depth;

    /**
     * Read counts is the kmer counts in each collection, the kmers are in the same order as in graph_kmers.
     */
    std::vector<std::vector<uint16_t>> read_counts;
    uint64_t max_mem;

    /**
     * Frequency value of the main unique content frequency peak for the first read_count collection.
     * This value usually is NOT manipulated by hand
     */
    uint16_t uniq_mode=0;

    static const bsgVersion_t min_compat;
};


#endif //BSG_KMERCOMPRESSIONINDEX_HPP
