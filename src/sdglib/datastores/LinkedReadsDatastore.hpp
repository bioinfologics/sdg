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
#include "ReadSequenceBuffer.hpp"

enum class LinkedReadsFormat {UCDavis,raw};
struct LinkedReadData {
    LinkedTag tag;
    std::string seq1,seq2;

    inline bool operator<(const struct LinkedReadData &other) const{
        if (tag<other.tag) return true;
        return false;
    }
};

std::string LinkedTag_to_seq(LinkedTag tag, uint8_t k=16);

/**
 * LinkedReadsDatastore is a file with reads and 10xTags.
 * Reads are stored sorted by tag, to optimise access when performing tag-based analyses
 * The binary file contains the following data:
 *
 * uint64_t readsize;
 * uint64_t read_tag.size()=pairs; (tag 0 is for pair (1,2) and 1 for (2,3), etc)
 * bsg10xTag[pairs] contents of rhe read_tag vector, a tag for each read pair.
 * read sequences: as \\0 terminated characters, using 2*readsize+2 for each pair
 *
 *
 */
class LinkedReadsDatastore {
public:
    LinkedReadsDatastore(WorkSpace &ws, std::string filename);
    LinkedReadsDatastore(WorkSpace &ws, std::string filename, std::ifstream &infile);
    LinkedReadsDatastore(WorkSpace &ws);
    LinkedReadsDatastore(WorkSpace &ws, std::string read1_filename, std::string read2_filename,
                         std::string output_filename, LinkedReadsFormat format, std::string default_name="",
                         int readsize = 250);
    LinkedReadsDatastore(WorkSpace &ws, const LinkedReadsDatastore &o);

    ~LinkedReadsDatastore(){
        if (fd) {
            fclose(fd);
        }
    }


    /**
     * @brief Provides an overview of the information in the LinkedReadsDatastore
     * @param level Base indentation level to use on the result
     * @param recursive Whether it should explore or not the rest of the hierarchy
     * @return
     * A text summary of the information contained in a LinkedReadsDatastore
    */
    std::string ls(int level=0,bool recursive=true) const;

    friend std::ostream& operator<<(std::ostream& os, const LinkedReadsDatastore &lrds);

    void print_status() const;
    static void build_from_fastq(std::string output_filename, std::string default_name, std::string read1_filename, std::string read2_filename,
                                 LinkedReadsFormat format, uint64_t readsize = 250, size_t chunksize = 10000000);
    void write(std::ofstream & output_file);
    void read(std::ifstream & input_file);
    void load_index(std::string _filename);
    //void read_index(std::ifstream & input_file);

    size_t size() const {return read_tag.size()*2;};
    std::string get_read_sequence(size_t readID);
    size_t get_read_pair(size_t readID) {return (readID%2==1) ? readID+1: readID-1;};
    LinkedTag get_read_tag(size_t readID);
    std::unordered_set<uint64_t> get_tags_kmers(int k, int min_tag_cov, std::set<LinkedTag> tags, ReadSequenceBuffer & blrsg);
    std::unordered_set<__uint128_t, int128_hash> get_tags_kmers128(int k, int min_tag_cov, std::set<LinkedTag> tags, ReadSequenceBuffer & blrsg, bool count_tag_cvg=false);
    std::vector<uint64_t> get_tag_reads(LinkedTag tag) const;
    std::vector<std::pair<LinkedTag, uint32_t>> get_tag_readcount();
    void dump_tag_occupancy_histogram(std::string filename);
    std::string filename; //if store is in single file sdg format these two are the same as the index file.

    uint64_t readsize;
    uint64_t readpos_offset;
    LinkedReadsMapper mapper;
    std::string name;
    std::string default_name;
private:
    std::vector<uint32_t> read_tag;
    FILE * fd=NULL;
    static const sdgVersion_t min_compat;
    WorkSpace &ws;


    //TODO: read sequence cache (std::map with a limit of elements and use count)
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
    std::unordered_set<uint64_t> get_tags_kmers(int min_tag_cov, std::set<LinkedTag> tags);
    void get_tag_kmers(LinkedTag tag);

private:

    struct tag_kmers_t {
        LinkedTag tag;
        std::unordered_set<uint64_t> kmers;
        const bool operator==(const LinkedTag & other_tag){return tag==other_tag;};
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
    ReadSequenceBuffer bprsg;
    StreamKmerFactory skf;

private:
    std::vector<uint64_t> counts;
};

