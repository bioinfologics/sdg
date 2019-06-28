//
// Created by Bernardo Clavijo (EI) on 2019-06-25.
//

#pragma once


#include <sdglib/workspace/WorkSpace.hpp>

class KmerCountsDatastore {
public:
    KmerCountsDatastore(const WorkSpace &_ws, uint8_t _k):ws(_ws),k(_k){
        index_sdg();
    };
    void index_sdg();


    /**
     * Accumulates the kmer count from the provided fastq file to the last available read_counts collection
     * @param filename Path to fastq file
     */
    //void add_count(std::vector<std::string> filenames);

    /**
     * Accumulates the kmer count from the provided data-store to the last available read_counts collection
     * @param ds PairedReadsDatastore ds
     */

    void add_count(const std::string & count_name, const PairedReadsDatastore & datastore);
    void add_count(const std::string & count_name, const LinkedReadsDatastore & datastore);
    void add_count(const std::string & count_name, const LongReadsDatastore & datastore);

    std::unordered_map<uint64_t,uint64_t> & get_kindex_map() const;

    std::vector<uint64_t> kindex;
    std::vector<std::string> count_names;
    std::vector<std::vector<uint16_t>> counts;
private:

    void add_string_to_count(int16_t count_idx,char * s);

    const WorkSpace &ws;
    const int8_t k;
};

