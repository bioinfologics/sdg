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

class BufferedPairedSequenceGetter;

class PairedReadsDatastore {
public:
    PairedReadsDatastore(WorkSpace &ws, std::ifstream &infile);
    PairedReadsDatastore(WorkSpace &ws, std::string filename, std::ifstream &infile);
    PairedReadsDatastore(WorkSpace &ws, std::string _filename);
    PairedReadsDatastore(WorkSpace &ws, std::string read1_filename,std::string read2_filename, std::string output_filename, int min_readsize=0, int max_readsize=250);
    PairedReadsDatastore(WorkSpace &ws, PairedReadsDatastore &ds);

    ~PairedReadsDatastore(){
        if (fd) {
            fclose(fd);
        }
    }

    PairedReadsDatastore& operator=(PairedReadsDatastore const &o);
    void print_status();
    static void build_from_fastq(std::string read1_filename,std::string read2_filename, std::string output_filename, int min_readsize=0, int max_readsize=250, size_t chunksize=10000000);
    void write(std::ofstream & output_file);
    void write_selection(std::ofstream &output_file, std::vector<uint64_t> read_ids);
    void read(std::ifstream & input_file);
    void load_index();
    uint64_t size()const;
    std::string get_read_sequence(size_t readID);
    std::unordered_set<__uint128_t, int128_hash> get_all_kmers128(int k, int min_tag_cov);
    std::unordered_set<__uint128_t, int128_hash> get_reads_kmers128(int k, int min_tag_cov, std::vector<uint64_t> reads);


    std::string filename; //if store is in single file sdg format these two are the same as the index file.
    uint64_t readsize;
    uint64_t readpos_offset;
    PairedReadsMapper mapper;
private:
    //TODO: save size
    uint64_t _size;
    FILE * fd=NULL;
    static const bsgVersion_t min_compat = 0x0001;

    WorkSpace &ws;

    //TODO: read sequence cache (std::map with a limit of elements and use count)
};

class BufferedPairedSequenceGetter{
public:
    BufferedPairedSequenceGetter(const PairedReadsDatastore &_ds, size_t _bufsize, size_t _chunk_size):
            datastore(_ds),bufsize(_bufsize),chunk_size(_chunk_size){
        fd=open(datastore.filename.c_str(),O_RDONLY);
        buffer=(char *)malloc(bufsize);
        buffer_offset=SIZE_MAX;
    }
    const char * get_read_sequence(uint64_t readID);
    ~BufferedPairedSequenceGetter(){
        free(buffer);
        close(fd);
    }
private:
    const PairedReadsDatastore &datastore;
    char * buffer;
    size_t bufsize,chunk_size;
    size_t buffer_offset;
    int fd;
};

