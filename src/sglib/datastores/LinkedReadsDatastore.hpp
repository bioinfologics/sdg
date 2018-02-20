//
// Created by Bernardo Clavijo (EI) on 10/02/2018.
//

#ifndef BSG_LINKEDREADSDATASTORE_HPP
#define BSG_LINKEDREADSDATASTORE_HPP

#include <sglib/PairedReadMapper.h>
#include <cstddef>

typedef uint32_t bsg10xTag;
enum LinkedReadsFormat {UCDavis,raw,seq};
unsigned const group_size=30;


class LinkedReadsDatastore {
public:
    LinkedReadsDatastore(std::string read1_filename,std::string read2_filename,LinkedReadsFormat format){
        build_from_fastq(read1_filename,read2_filename,format);
    };
    void build_from_fastq(std::string read1_filename,std::string read2_filename,LinkedReadsFormat format);
    size_t size(){return read_offset.size()-1;};

    std::string get_read_sequence(size_t readID,FILE * file1, FILE * file2);
    std::string get_read_sequence_fd(size_t readID,int fd1, int fd2);
    void get_read_sequence_fd(size_t readID, int fd1, int fd2, char * dest);
    inline std::string get_read_sequence(size_t readID){return get_read_sequence(readID,fd1,fd2);};
    bsg10xTag get_read_tag(size_t readID);
    void dump_index_to_disk(std::string filename);
    void dump_full_store_to_disk(std::string filename);
    void load_from_disk(std::string filename,std::string read1_filename="",std::string read2_filename="");
    std::string filename1,filename2; //if store is in single file bsg format these two are the same as the index file.

    std::vector<uint64_t> group_offset1;
    std::vector<uint64_t> group_offset2;
    std::vector<uint16_t> read_offset;
private:
    std::vector<uint32_t> read_tag;
    FILE * fd1=NULL;
    FILE * fd2=NULL;
    //TODO: read sequence cache (std::map with a limit of elements and use count)
};

class BufferedLRSequenceGetter{
public:
    BufferedLRSequenceGetter(const LinkedReadsDatastore &_ds, size_t _bufsize, size_t _chunk_size):
            datastore(_ds),bufsize(_bufsize),chunk_size(_chunk_size){
        buffer1=(char *)malloc(bufsize);
        buffer2=(char *)malloc(bufsize);
        fd1=open(datastore.filename1.c_str(),O_RDONLY);
        fd2=open(datastore.filename2.c_str(),O_RDONLY);
        buffer1_offset=SIZE_MAX;
        buffer1_offset=SIZE_MAX;
    }
    const char * get_read_sequence(uint64_t readID);
    ~BufferedLRSequenceGetter(){
        free(buffer1);
        free(buffer2);
        close(fd1);
        close(fd2);
    }
private:
    const LinkedReadsDatastore &datastore;
    char * buffer1;
    char * buffer2;
    size_t bufsize,chunk_size;
    size_t buffer1_offset,buffer2_offset;
    int fd1,fd2;
};

#endif //BSG_LINKEDREADSDATASTORE_HPP

