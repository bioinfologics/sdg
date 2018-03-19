//
// Created by Bernardo Clavijo (EI) on 10/02/2018.
//

#ifndef BSG_LINKEDREADSDATASTORE_HPP
#define BSG_LINKEDREADSDATASTORE_HPP

#include <sglib/PairedReadMapper.hpp>
#include <cstddef>

typedef uint32_t bsg10xTag;
enum LinkedReadsFormat {UCDavis,raw,seq};
struct readData {
    bsg10xTag tag;
    std::string seq1,seq2;

    inline bool operator<(const struct readData &other) const{
        if (tag<other.tag) return true;
        return false;
    }
};
class BufferedLRSequenceGetter;

class LinkedReadsDatastore {
public:
    LinkedReadsDatastore(){};
    LinkedReadsDatastore(std::string filename){
        load_index(filename);
    };
    LinkedReadsDatastore(std::string read1_filename,std::string read2_filename, std::string output_filename, LinkedReadsFormat format, int readsize=250){
        build_from_fastq(read1_filename,read2_filename,output_filename,format,readsize);
    };
    void build_from_fastq(std::string read1_filename,std::string read2_filename, std::string output_filename, LinkedReadsFormat format, int readsize=250,size_t chunksize=10000000);
    void write(std::ofstream & output_file);
    void read(std::ifstream & input_file);
    void load_index(std::string _filename);
    //void read_index(std::ifstream & input_file);

    size_t size(){return read_tag.size()*2-1;};
    std::string get_read_sequence(size_t readID);
    //inline std::string get_read_sequence(size_t readID){return get_read_sequence(readID,fd1,fd2);};
    bsg10xTag get_read_tag(size_t readID);
    std::unordered_set<uint64_t> get_tags_kmers(int k, int min_tag_cov, std::unordered_set<bsg10xTag> tags, BufferedLRSequenceGetter & blrsg);
    std::vector<uint64_t> get_tag_reads(bsg10xTag tag);
    void dump_tag_occupancy_histogram(std::string filename);
    std::string filename; //if store is in single file bsg format these two are the same as the index file.

    uint64_t readsize;
    uint64_t readpos_offset;
private:
    std::vector<uint32_t> read_tag;
    FILE * fd=NULL;

    //TODO: read sequence cache (std::map with a limit of elements and use count)
};

class BufferedLRSequenceGetter{
public:
    BufferedLRSequenceGetter(const LinkedReadsDatastore &_ds, size_t _bufsize, size_t _chunk_size):
            datastore(_ds),bufsize(_bufsize),chunk_size(_chunk_size){
        fd=open(datastore.filename.c_str(),O_RDONLY);
        buffer=(char *)malloc(bufsize);
        buffer_offset=SIZE_MAX;
    }
    const char * get_read_sequence(uint64_t readID);
    ~BufferedLRSequenceGetter(){
        free(buffer);
        close(fd);
    }
private:
    const LinkedReadsDatastore &datastore;
    char * buffer;
    size_t bufsize,chunk_size;
    size_t buffer_offset;
    int fd;
};

#endif //BSG_LINKEDREADSDATASTORE_HPP

