//
// Created by Bernardo Clavijo (EI) on 2019-06-25.
//

#pragma once

#include <string>
#include <vector>
#include <tuple>
#include <fstream>
#include <unordered_map>

enum KmerCountMode{Canonical,NonCanonical};

class WorkSpace;
class PairedReadsDatastore;
class LinkedReadsDatastore;
class LongReadsDatastore;

/**
 * KmerCounters are containers for the coverage of a set of kmers from a graph or workspace,
 * they provide functionality to count and store an index of the kmers present and
 * multiple instances of counts for these kmers.
 */
class KmerCounter {
public:
    KmerCounter (const WorkSpace &_ws, const std::string &_name, uint8_t _k, KmerCountMode _count_mode=Canonical):ws(_ws),k(_k), name(_name),count_mode(_count_mode){
        index_sdg();
    };
    KmerCounter (const WorkSpace &ws, std::ifstream &infile);
    KmerCounter (const WorkSpace &_ws, const std::string &filename):ws(_ws) {
        std::ifstream count_file(filename);
        read_counts(count_file);
    }


    /**
     * @brief Provides an overview of the information in the KmerCounter
     * @param level Base indentation level to use on the result
     * @param recursive Whether it should explore or not the rest of the hierarchy
     * @return
     * A text summary of the information contained in a KmerCounter
    */
    std::string ls(int level=0,bool recursive=true) const;

    friend std::ostream& operator<<(std::ostream &os, const KmerCounter &kc);

    void index_sdg();

    void update_graph_counts();

    bool operator==(const KmerCounter &o) const {
        return (std::tie(k, kindex, count_names, counts) == std::tie(o.k, o.kindex, o.count_names, o.counts));
    }

    KmerCounter& operator=(const KmerCounter &o) {
        if ( &o == this) return *this;

        counts = o.counts;
        count_names = o.count_names;
        index_sdg();
    }

    /**
     * @brief Accumulates the kmer count from the provided fastq file to the last available read_counts collection
     * @param filename Path to fastq file
     */
    void add_count(const std::string & count_name,const std::vector<std::string> &filenames, bool fastq=true);

    /**
     * @brief Accumulates the kmer count from the provided data-store to the last available read_counts collection
     * @param ds PairedReadsDatastore ds
     */
    void add_count(const std::string & count_name, const PairedReadsDatastore & datastore);

    /**
     * @brief Accumulates the kmer count from the provided data-store to the last available read_counts collection
     * @param ds LinkedReadsDatastore ds
     */
    void add_count(const std::string & count_name, const LinkedReadsDatastore & datastore);

    /**
     * @brief Accumulates the kmer count from the provided data-store to the last available read_counts collection
     * @param ds LongReadsDatastore ds
     */
    void add_count(const std::string & count_name, const LongReadsDatastore & datastore);

    /**TODO
     * @brief Accumulates the kmer count from the DG's nodes, allowing to filter by size and connection status
     * @param ds LongReadsDatastore ds
     */
    //void add_count(const std::string & count_name, const DistanceGraph & dg, uint32_t min_size=0, bool only_connected=true);

    /**
     * @brief Retrieves the counts for each kmer in sequence s from count_name
     * @param count_name Name of the count to query
     * @param s Sequence to project
     * @return Vector of kmer counts for each kmer in s
     */
    std::vector<uint16_t> project_count(const std::string & count_name, const std::string &s);

    /**
     * @brief Retrieves the counts for each kmer in sequence s from counts[count_idx]
     * @param count_idx Index of the count to query
     * @param s Sequence to project
     * @return Vector of kmer counts for each kmer in s
     */
    std::vector<uint16_t> project_count(const uint16_t count_idx, const std::string &s);

    void set_kci_peak(float f){
        if (f!=kci_peak_f){
            kci_peak_f=f;
            kci_cache.clear();
        }

    }

    float kci(int64_t node);

    std::vector<uint64_t> count_spectra(std::string name, uint16_t maxf=1000, bool unique_in_graph=false);

    void write(std::ofstream & output_file) const;
    void write(std::fstream & output_file) const;
    void write_counts(std::ofstream &count_file) const;

    void read(std::ifstream & input_file);
    void read_counts(std::ifstream &count_file);
    int8_t get_k(){return k;};


    /**
     * @brief Retrieves a count by name
     * @param name Name of the count
     * @return Vector of counts for each kmer in the KmerCounter index
     */
    const std::vector<uint16_t> & get_count_by_name(const std::string &name) const;
    std::vector<std::string> list_names ();

    std::vector<uint64_t> kindex;               /// Ordered list of kmers that contain counts
    std::vector<std::string> count_names;       /// Names of the counts vectors
    std::vector<std::vector<uint16_t>> counts;  /// Count vector, contains an entry per kmer in the kindex
    
    std::string name;   /// Name of the KmerCounter

private:
    const WorkSpace &ws;
    int8_t k;
    KmerCountMode count_mode;
    float kci_peak_f=-1;
    std::unordered_map<int64_t, float> kci_cache;
};

