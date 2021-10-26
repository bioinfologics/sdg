//
// Created by Bernardo Clavijo (EI) on 10/05/2018.
//

#pragma once

#include <memory>
#include <cstdint>
#include <iostream>
#include <string>
#include <unordered_set>
#include <set>
#include <vector>
#include <cstddef>
#include <list>
#include <fcntl.h>
#include <unistd.h>
#include <algorithm>
#include <sdglib/factories/KMerFactory.hpp>
#include <sdglib/mappers/PairedReadsMapper.hpp>
#include "sdglib/Version.hpp"

class WorkSpace;
struct PairedReadData {
    std::string seq1,seq2;
};

class PairedReadsDatastore {
public:
    PairedReadsDatastore(WorkSpace &ws);
    PairedReadsDatastore(WorkSpace &ws, std::string filename, std::ifstream &infile);
    PairedReadsDatastore(WorkSpace &ws, std::string _filename);
    PairedReadsDatastore(WorkSpace &ws, std::string read1_filename,std::string read2_filename, std::string output_filename, std::string default_name="", int min_readsize=0, int max_readsize=250, int fs=0, int orientation=0);

    ~PairedReadsDatastore(){
        if (fd) {
            fclose(fd);
        }
    }

    /**
     * @brief Provides an overview of the information in the PairedReadsDatastore
     * @param level Base indentation level to use on the result
     * @param recursive Whether it should explore or not the rest of the hierarchy
     * @return
     * A text summary of the information contained in a PairedReadsDatastore
    */
    std::string ls(int level=0,bool recursive=true) const;

    friend std::ostream& operator<<(std::ostream &os, const PairedReadsDatastore &prds);
    void print_status() const;
    static void build_from_fastq(std::string output_filename, std::string read1_filename, std::string read2_filename, std::string default_name,
                                 uint64_t min_readsize = 0, uint64_t max_readsize = 250, int fragment_size=0, int orientation=0, size_t chunksize = 10000000, int run_length_limit=0);
    void write(std::ofstream & output_file);
    void read(std::ifstream & input_file);
    void load_index();
    uint64_t size()const;
    std::string get_read_sequence(seqID_t readID);
    seqID_t get_read_pair(seqID_t readID) const;
    static uint32_t read_sanitizer(char* bp);
    std::unordered_set<__uint128_t, int128_hash> get_all_kmers128(int k, int min_tag_cov);
    std::unordered_set<__uint128_t, int128_hash> get_reads_kmers128(int k, int min_tag_cov, std::vector<uint64_t> reads);


    std::string filename; //if store is in single file sdg format these two are the same as the index file.
    uint64_t readsize;
    uint64_t readpos_offset;

    uint64_t fragment_size;
    uint64_t orientation; // 0: Undefined, 1: FW-REV (PE), 2: REV-FW (MP)

    PairedReadsMapper mapper;
    std::string name;
    std::string default_name;
private:
    //TODO: save size
    uint64_t _size;
    FILE * fd=NULL;
    static const sdgVersion_t min_compat = 0x0003;

    WorkSpace &ws;

    //TODO: read sequence cache (std::map with a limit of elements and use count)
};

