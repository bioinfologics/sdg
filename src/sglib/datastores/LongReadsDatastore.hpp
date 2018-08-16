//
// Created by Luis Yanes (EI) on 23/03/2018.
//

#ifndef BSG_LONGREADSDATASTORE_HPP
#define BSG_LONGREADSDATASTORE_HPP

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
#include <sglib/Version.hpp>

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
    int lr_sequence_fd = 0;
    std::string file_containing_long_read_sequence;

    void read_rtfr(std::ifstream &ifs);
    void write_rtfr(std::ofstream &output_file);
    void load_index(std::string &file);


public:
    std::vector< ReadPosSize > read_to_fileRecord{ReadPosSize(0,0)};
    /**
     * Initialize from already created index
     * @param filename
     *
     * Initialises the memory mapping of the reads file
     */
    explicit LongReadsDatastore(std::string filename);
    LongReadsDatastore() = default;
    /**
     * Initialize from long_read_file then store the index
     * @param long_read_file
     * @param output_file
     *
     * Initialises the memory mapping of the reads file
     */
    LongReadsDatastore(std::string long_read_file, std::string output_file);
    uint32_t build_from_fastq(std::ofstream &outf, std::string long_read_file);

    void read(std::ifstream &ifs);
    void write(std::ofstream &output_file);

    size_t size() const { return read_to_fileRecord.size(); }

    void load_from_stream(std::string file_name, std::ifstream &input_file);
    std::string filename;
    static const bsgVersion_t min_compat;
};

// Check if this needs to be page size aware
class BufferedSequenceGetter{
public:
    BufferedSequenceGetter(const LongReadsDatastore &_ds, size_t _bufsize = 1024*1024*30, size_t _chunk_size = 1024*1024*4):
            datastore(_ds),bufsize(_bufsize),chunk_size(_chunk_size){
        fd=open(datastore.filename.c_str(),O_RDONLY);
        if (fd == -1) {
            std::string msg("Cannot open file " + datastore.filename);
            perror(msg.c_str());
            throw std::runtime_error("Cannot open " + datastore.filename);
        }

        buffer=(char *)malloc(bufsize);
    }
    std::string get_read_sequence(uint64_t readID);

    ~BufferedSequenceGetter(){
        free(buffer);
        close(fd);
    }

    void write_selection(std::ofstream &output_file, const std::vector<uint64_t> &read_ids);

private:
    const LongReadsDatastore &datastore;
    char * buffer;
    size_t bufsize,chunk_size;
    off_t buffer_offset = std::numeric_limits<off_t>::max();
    int fd;
};

#endif //BSG_LONGREADSDATASTORE_HPP
