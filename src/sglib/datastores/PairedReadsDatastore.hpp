//
// Created by Bernardo Clavijo (EI) on 10/05/2018.
//

#ifndef BSG_PAIREDEDREADSDATASTORE_HPP
#define BSG_PAIREDEDREADSDATASTORE_HPP

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
#include <sglib/factories/KMerFactory.hpp>
#include "sglib/Version.hpp"


struct PairedReadData {
    std::string seq1,seq2;
};

class BufferedPairedSequenceGetter;

class PairedReadsDatastore {
public:
    PairedReadsDatastore(){};
    PairedReadsDatastore(std::string _filename){
        filename=_filename;
        load_index();
    };
    PairedReadsDatastore(std::string read1_filename,std::string read2_filename, std::string output_filename, int min_readsize=0, int max_readsize=250){
        build_from_fastq(read1_filename,read2_filename,output_filename,min_readsize,max_readsize);
    };
    void build_from_fastq(std::string read1_filename,std::string read2_filename, std::string output_filename, int min_readsize=0, int max_readsize=250, size_t chunksize=10000000);
    void write(std::ofstream & output_file);
    void write_selection(std::ofstream &output_file, std::vector<uint64_t> read_ids);
    void read(std::ifstream & input_file);
    void load_index();
    void load_from_stream(std::string filename,std::ifstream & input_file);
    uint64_t size()const {return _size;};
    std::string get_read_sequence(size_t readID);
    std::unordered_set<__uint128_t, int128_hash> get_all_kmers128(int k, int min_tag_cov);
    std::unordered_set<__uint128_t, int128_hash> get_reads_kmers128(int k, int min_tag_cov, std::vector<uint64_t> reads);
    std::string filename; //if store is in single file bsg format these two are the same as the index file.

    uint64_t readsize;
    uint64_t readpos_offset;
private:
    //TODO: save size
    uint64_t _size;
    FILE * fd=NULL;
    static const bsgVersion_t min_compat = 0x0001;

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

#endif //BSG_LINKEDREADSDATASTORE_HPP

