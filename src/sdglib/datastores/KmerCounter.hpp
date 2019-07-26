//
// Created by Bernardo Clavijo (EI) on 2019-06-25.
//

#pragma once

#include <string>
#include <vector>
#include <tuple>
#include <fstream>

enum KmerCountMode{Canonical,NonCanonical};

class WorkSpace;
class PairedReadsDatastore;
class LinkedReadsDatastore;
class LongReadsDatastore;
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

    std::string ls(int level=0,bool recursive=true);

    void index_sdg();

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
     * Accumulates the kmer count from the provided fastq file to the last available read_counts collection
     * @param filename Path to fastq file
     */
    void add_count(const std::string & count_name,const std::vector<std::string> &filenames, bool fastq=true);

    /**
     * Accumulates the kmer count from the provided data-store to the last available read_counts collection
     * @param ds PairedReadsDatastore ds
     */
    void add_count(const std::string & count_name, const PairedReadsDatastore & datastore);
    void add_count(const std::string & count_name, const LinkedReadsDatastore & datastore);
    void add_count(const std::string & count_name, const LongReadsDatastore & datastore);

    std::vector<uint16_t> project_count(const std::string & count_name, const std::string &s);
    std::vector<uint16_t> project_count(const uint16_t count_idx, const std::string &s);

    void write(std::ofstream & output_file) const;
    void write(std::fstream & output_file) const;
    void write_counts(std::ofstream &count_file) const;

    void read(std::ifstream & input_file);
    void read_counts(std::ifstream &count_file);
    int8_t get_k(){return k;};


    const std::vector<uint16_t> & get_count_by_name(const std::string &name) const;
    std::vector<std::string> get_count_names ();

    std::vector<uint64_t> kindex;
    std::vector<std::string> count_names;
    std::vector<std::vector<uint16_t>> counts;
    
    std::string name;

private:
    const WorkSpace &ws;
    int8_t k;
    KmerCountMode count_mode;
};

