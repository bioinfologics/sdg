//
// Created by Bernardo Clavijo (EI) on 10/02/2018.
//

#ifndef BSG_LINKEDREADSDATASTORE_HPP
#define BSG_LINKEDREADSDATASTORE_HPP

#include <sglib/PairedReadMapper.hpp>
#include <cstddef>

typedef uint32_t bsg10xTag;
enum LinkedReadsFormat {UCDavis,raw};

class LinkedRead{
    size_t offset;
    bsg10xTag tag;
    uint16_t size;
};

class LinkedReadsDatastore {
public:
    void build_from_fastq(std::string read1_filename,std::string read2_filename,LinkedReadsFormat format);
    std::string get_read_sequence(size_t readID);
    void dump_index_to_disk(std::string filename);
    void dump_full_store_to_disk(std::string filename);
    void load_from_disk(std::string filename,std::string read1_filename="",std::string read2_filename="");
private:
    std::string filename1,filename2; //if store is in single file bsg format these two are the same as the index file.
    std::vector<LinkedRead> reads;
    std::unordered_map<bsg10xTag,std::vector<size_t>> reads_in_tag;
    //TODO: read sequence cache (std::map with a limit of elements and use count)
};


#endif //BSG_LINKEDREADSDATASTORE_HPP

