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

    std::string file_containing_long_read_sequence;

    void load_index(std::string &file);

    WorkSpace &ws;

public:
    std::vector< ReadPosSize > read_to_fileRecord{ReadPosSize(0,0)};

    LongReadsDatastore(WorkSpace &ws, std::ifstream &infile);
    LongReadsDatastore(WorkSpace &ws, std::string filename, std::ifstream &input_file);
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
    LongReadsDatastore(WorkSpace &ws, std::string long_read_file, std::string output_file);

    LongReadsDatastore& operator=(LongReadsDatastore const &o);
    uint32_t build_from_fastq(std::ofstream &outf, std::string long_read_file);
    static void build_from_fastq(std::string outf, std::string long_read_file);
    void print_status();
    void read(std::ifstream &ifs);
    void write(std::ofstream &output_file);

    size_t size() const { return read_to_fileRecord.size(); }

    std::string filename;
    static const bsgVersion_t min_compat;

    LongReadsMapper mapper;
};

// Check if this needs to be page size aware
class BufferedSequenceGetter{
public:
    BufferedSequenceGetter(const LongReadsDatastore &_ds, size_t _bufsize = (1024*1024*30ul), size_t _chunk_size = (1024*1024*4ul)):
            datastore(_ds),bufsize(_bufsize),chunk_size(_chunk_size){
        fd=open(datastore.filename.c_str(),O_RDONLY);
        if (fd == -1) {
            std::string msg("Cannot open file " + datastore.filename);
            perror(msg.c_str());
            throw std::runtime_error("Cannot open " + datastore.filename);
        }
        struct stat f_stat;
        stat(_ds.filename.c_str(), &f_stat);
        total_size = f_stat.st_size;
        buffer=(char *)malloc(bufsize);
    }
    const char * get_read_sequence(uint64_t readID);

    ~BufferedSequenceGetter(){
        free(buffer);
        if(fd) close(fd);
    }

    void write_selection(std::ofstream &output_file, const std::vector<uint64_t> &read_ids);

private:
    const LongReadsDatastore &datastore;
    char * buffer;
    size_t bufsize,chunk_size;
    off_t buffer_offset = std::numeric_limits<off_t>::max();
    int fd;
    off_t total_size=0;
};
