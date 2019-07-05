//
// Created by Bernardo Clavijo (EI) on 2019-06-25.
//

#pragma once

#include <string>
#include <vector>

class WorkSpace;
class PairedReadsDatastore;
class LinkedReadsDatastore;
class LongReadsDatastore;
class KmerCountsDatastore {
public:
    KmerCountsDatastore(const WorkSpace &_ws, const std::string &_name, uint8_t _k):ws(_ws),k(_k), name(_name){
        index_sdg();
    };
    KmerCountsDatastore (const WorkSpace &ws, std::ifstream &infile);
    void index_sdg();

    KmerCountsDatastore& operator=(const KmerCountsDatastore &o) {
        if ( &o == this) return *this;

        counts = o.counts;
        count_names = o.count_names;
        index_sdg();
    }

    /**
     * Accumulates the kmer count from the provided fastq file to the last available read_counts collection
     * @param filename Path to fastq file
     */
    void add_count(const std::string & count_name,const std::vector<std::string> &filenames);

    /**
     * Accumulates the kmer count from the provided data-store to the last available read_counts collection
     * @param ds PairedReadsDatastore ds
     */
    void add_count(const std::string & count_name, const PairedReadsDatastore & datastore);
    void add_count(const std::string & count_name, const LinkedReadsDatastore & datastore);
    void add_count(const std::string & count_name, const LongReadsDatastore & datastore);

    std::vector<uint16_t> project_count(const std::string & count_name, const std::string &s);
    std::vector<uint16_t> project_count(const uint16_t count_idx, const std::string &s);

    void write(std::ofstream & output_file);
    void read(std::ifstream & input_file);
    int8_t get_k(){return k;};

    std::vector<uint64_t> kindex;
    std::vector<std::string> count_names;
    std::vector<std::vector<uint16_t>> counts;
    
    std::string name;

private:
    const WorkSpace &ws;
    int8_t k;
};

