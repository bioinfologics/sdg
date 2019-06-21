//
// Created by Bernardo Clavijo (EI) on 10/02/2018.
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
#include <sdglib/Version.hpp>
#include <sdglib/mappers/LinkedReadsMapper.hpp>

enum class LinkedReadsFormat {UCDavis,raw,seq};
struct LinkedReadData {
    bsg10xTag tag;
    std::string seq1,seq2;

    inline bool operator<(const struct LinkedReadData &other) const{
        if (tag<other.tag) return true;
        return false;
    }
};
class BufferedLRSequenceGetter;

std::string bsg10xTag_to_seq(bsg10xTag tag, uint8_t k=16);

/**
 * LinkedReadsDatastore is a file with reads and 10xTags.
 * Reads are stored sorted by tag, to optimise access when performing tag-based analyses
 * The binary file contains the following data:
 *
 * uint64_t readsize;
 * uint64_t read_tag.size()=pairs; (tag 0 is for pair (1,2) and 1 for (2,3), etc)
 * bsg10xTag[pairs] contents of rhe read_tag vector, a tag for each read pair.
 * read sequences: as \0 terminated characters, using 2*readsize+2 for each pair
 *
 *
 */
class LinkedReadsDatastore {
public:
    LinkedReadsDatastore(WorkSpace &ws,std::string filename);
    LinkedReadsDatastore(WorkSpace &ws,std::string filename, std::ifstream &infile);
    LinkedReadsDatastore(WorkSpace &ws, std::ifstream &infile);
    LinkedReadsDatastore(WorkSpace &ws,std::string read1_filename,std::string read2_filename, std::string output_filename, LinkedReadsFormat format, int readsize=250);
    LinkedReadsDatastore(WorkSpace &ws, const LinkedReadsDatastore &o);

    ~LinkedReadsDatastore(){
        if (fd) {
            fclose(fd);
        }
    }

    LinkedReadsDatastore& operator=(LinkedReadsDatastore const &o);

    void print_status();
    static void build_from_fastq(std::string output_filename, std::string read1_filename, std::string read2_filename,
                                 LinkedReadsFormat format, uint64_t readsize = 250, size_t chunksize = 10000000);
    void write(std::ofstream & output_file);
    void write_selection(std::ofstream & output_file, const std::set<bsg10xTag> & tagSet);
    void read(std::ifstream & input_file);
    void load_index(std::string _filename);
    //void read_index(std::ifstream & input_file);

    size_t size() const {return read_tag.size()*2;};
    std::string get_read_sequence(size_t readID);
    //inline std::string get_read_sequence(size_t readID){return get_read_sequence(readID,fd1,fd2);};
    bsg10xTag get_read_tag(size_t readID);
    std::unordered_set<uint64_t> get_tags_kmers(int k, int min_tag_cov, std::set<bsg10xTag> tags, BufferedLRSequenceGetter & blrsg);
    std::unordered_set<__uint128_t, int128_hash> get_tags_kmers128(int k, int min_tag_cov, std::set<bsg10xTag> tags, BufferedLRSequenceGetter & blrsg, bool count_tag_cvg=false);
    std::vector<uint64_t> get_tag_reads(bsg10xTag tag) const;
    std::vector<std::pair<bsg10xTag, uint32_t>> get_tag_readcount();
    void dump_tag_occupancy_histogram(std::string filename);
    std::string filename; //if store is in single file sdg format these two are the same as the index file.

    uint64_t readsize;
    uint64_t readpos_offset;
    LinkedReadsMapper mapper;
private:
    std::vector<uint32_t> read_tag;
    FILE * fd=NULL;
    static const bsgVersion_t min_compat;
    WorkSpace &ws;


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
        if(fd) close(fd);
    }
private:
    const LinkedReadsDatastore &datastore;
    char * buffer;
    size_t bufsize,chunk_size;
    size_t buffer_offset;
    int fd;
};

/**
 * @brief kmerises sets of tags, but saves the kmers of individual tags in a buffer to speed-up going through many nodes
 *
 * @todo add an option to only use tags with X+ reads.
 */
class BufferedTagKmerizer{

public:
    BufferedTagKmerizer(const LinkedReadsDatastore &_ds, char K, size_t _bufsize, size_t _chunk_size):
            K(K),bprsg(_ds,_bufsize,_chunk_size),skf(K),datastore(_ds){counts.reserve(1000000);};
    std::unordered_set<uint64_t> get_tags_kmers(int min_tag_cov, std::set<bsg10xTag> tags);
    void get_tag_kmers(bsg10xTag tag);

private:

    struct tag_kmers_t {
        bsg10xTag tag;
        std::unordered_set<uint64_t> kmers;
        const bool operator==(const bsg10xTag & other_tag){return tag==other_tag;};
    };

    class StreamKmerFactory : public  KMerFactory {
    public:
        explicit StreamKmerFactory(uint8_t k) : KMerFactory(k){}
        inline void produce_all_kmers(const char * seq, std::vector<uint64_t> &mers){
            // TODO: Adjust for when K is larger than what fits in uint64_t!
            last_unknown=0;
            fkmer=0;
            rkmer=0;
            auto s=seq;
            while (*s!='\0' and *s!='\n') {
                //fkmer: grows from the right (LSB)
                //rkmer: grows from the left (MSB)
                fillKBuf(*s, fkmer, rkmer, last_unknown);
                if (last_unknown >= K) {
                    if (fkmer <= rkmer) {
                        // Is fwd
                        mers.emplace_back(fkmer);
                    } else {
                        // Is bwd
                        mers.emplace_back(rkmer);
                    }
                }
                ++s;
            }
        }
    };

    const LinkedReadsDatastore & datastore;
    uint8_t K;
    BufferedLRSequenceGetter bprsg;
    StreamKmerFactory skf;

private:
    std::vector<uint64_t> counts;
};

