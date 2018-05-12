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
#include <sglib/factories/KMerFactory.h>


struct PairedReadData {
    std::string seq1,seq2;
};

class BufferedPairedSequenceGetter;

class PairedReadsDatastore {
public:
    PairedReadsDatastore(){};
    PairedReadsDatastore(std::string filename){
        load_index(filename);
    };
    PairedReadsDatastore(std::string read1_filename,std::string read2_filename, std::string output_filename, int readsize=250){
        build_from_fastq(read1_filename,read2_filename,output_filename,readsize);
    };
    void build_from_fastq(std::string read1_filename,std::string read2_filename, std::string output_filename, int readsize=250,size_t chunksize=10000000);
    void write(std::ofstream & output_file);
    void read(std::ifstream & input_file);
    void load_index(std::string _filename);
    size_t size(){return _size;};
    std::string get_read_sequence(size_t readID);
    std::string filename; //if store is in single file bsg format these two are the same as the index file.

    uint64_t readsize;
    uint64_t readpos_offset;
private:
    //TODO: save size
    size_t _size;
    FILE * fd=NULL;

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

