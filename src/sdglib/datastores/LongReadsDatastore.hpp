//
// Created by Luis Yanes (EI) on 23/03/2018.
//

#pragma once

#include <memory>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>
#include <sys/mman.h>

#include <string>
#include <vector>
#include <fstream>
#include <unordered_map>
#include <cstring>

#include <sys/stat.h>
#include <limits>
#include <sdglib/Version.hpp>
#include <sdglib/mappers/LongReadsMapper.hpp>
#include "ReadSequenceBuffer.hpp"

class WorkSpace;
class LongReadsMapper;

struct ReadPosSize {
    off_t offset = 0;
    uint32_t record_size = 0;

    ReadPosSize() = default;
    ReadPosSize(off_t offset, uint32_t record_size) : offset(offset), record_size(record_size) {}

    bool operator==(const ReadPosSize &other) const {
        return (std::tie(offset, record_size) == std::tie(other.offset, other.record_size));
    }
};

class LongReadsDatastore {

    void load_index(std::string &file);

    WorkSpace &ws;
    int fd=NULL;
public:
    std::vector< ReadPosSize > read_to_fileRecord{ReadPosSize(0,0)};

    ~LongReadsDatastore(){
        if (fd) {
            close(fd);
        }
    }
    explicit LongReadsDatastore(WorkSpace &ws);
    LongReadsDatastore(WorkSpace &ws, std::string default_name, const std::string &filename, std::ifstream &input_file);

    /**
     * @brief Provides an overview of the information in the LongReadsDatastore
     * @param level Base indentation level to use on the result
     * @param recursive Whether it should explore or not the rest of the hierarchy
     * @return
     * A text summary of the information contained in a LongReadsDatastore
    */
    std::string ls(int level=0,bool recursive=true) const;

    /**
     * Initialize from already created index
     * @param filename
     *
     * Initialises the memory mapping of the reads file
     */
    LongReadsDatastore(WorkSpace &ws, std::string filename);

    friend std::ostream& operator<<(std::ostream &os, const LongReadsDatastore &lords);

    /**
     * Create a long reads data-store from fastq files and write the lords to disk
     * @param output_file Output filename of the datastore
     * @param default_name internal name of the datastore
     * @param long_read_file fastq file of long reads
     */
    static void build_from_fastq(const std::string &output_file, const std::string &default_name, const std::string &long_read_file, size_t min_size=0, int rle=0);
    void print_status() const;
    void read(std::ifstream &ifs);
    void write(std::ofstream &output_file);
    size_t size() const { return read_to_fileRecord.size()-1; }

    std::string get_read_sequence(size_t readID) const;

    std::string filename;
    std::string name;
    std::string default_name;
    static const sdgVersion_t min_compat;

    LongReadsMapper mapper;
};
