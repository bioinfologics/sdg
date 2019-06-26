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
    int fd;
public:
    std::vector< ReadPosSize > read_to_fileRecord{ReadPosSize(0,0)};

    ~LongReadsDatastore();
    LongReadsDatastore(WorkSpace &ws, std::ifstream &infile);
    LongReadsDatastore(WorkSpace &ws, const std::string &filename, std::ifstream &input_file);
    LongReadsDatastore(WorkSpace &ws, LongReadsDatastore &o);
    LongReadsDatastore(const LongReadsDatastore &o);
    /**
     * Initialize from already created index
     * @param filename
     *
     * Initialises the memory mapping of the reads file
     */
    LongReadsDatastore(WorkSpace &ws, std::string filename);
    /**
     * Initialize from long_read_file then store the index
     * @param long_read_file
     * @param output_file
     *
     * Initialises the memory mapping of the reads file
     */
    LongReadsDatastore(WorkSpace &ws, const std::string &long_read_file, const std::string &output_file);

    LongReadsDatastore& operator=(LongReadsDatastore const &o);
    uint32_t build_from_fastq(std::ofstream &outf, const std::string &long_read_file);
    static void build_from_fastq(const std::string &output_file, const std::string &long_read_file);
    void print_status() const;
    void read(std::ifstream &ifs);
    void write(std::ofstream &output_file);
    void write_selection(std::ofstream &output_file, const std::vector<uint64_t> &read_ids);
    size_t size() const { return read_to_fileRecord.size(); }

    std::string get_read_sequence(size_t readID) const;

    std::string filename;
    std::string name;
    static const sdgVersion_t min_compat;

    LongReadsMapper mapper;
};
