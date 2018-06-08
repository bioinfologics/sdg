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

#include <sglib/types/GenericTypes.hpp>
#include <sys/stat.h>

class LongReadsDatastore {
    struct ReadPosSize {
        off_t offset = 0;
        uint32_t record_size = 0;

        ReadPosSize() {}
        ReadPosSize(off_t offset, uint32_t record_size) : offset(offset), record_size(record_size) {}

        bool operator==(const ReadPosSize &other) const {
            return (std::tie(offset, record_size) == std::tie(other.offset, other.record_size));
        }
    };

    void *mapped_readfile = 0;
    off_t offset = 0;
    size_t file_size = 0;
    off_t page_size = static_cast<off_t>(sysconf(_SC_PAGESIZE));
    int long_read_fd = 0;
    std::string long_read_file;

    void read_rtfr(std::ifstream &ifs);
    void write_rtfr(std::ofstream &output_file);
    void load_index(std::string &file);


public:
    std::vector< ReadPosSize > read_to_fileRecord{ReadPosSize(0,0)};
    LongReadsDatastore(){};
    /**
     * Initialize from already created index
     * @param filename
     *
     * Initialises the memory mapping of the reads file
     */
    LongReadsDatastore(std::string filename);

    /**
     * Initialize from long_read_file then store the index
     * @param long_read_file
     * @param output_file
     *
     * Initialises the memory mapping of the reads file
     */
    LongReadsDatastore(std::string long_read_file, std::string output_file);
    uint32_t build_from_fastq(std::string long_read_file);

    void read(std::ifstream &ifs);
    void write(std::ofstream &output_file);

    size_t size() const { return read_to_fileRecord.size(); }
    std::string get_read_sequence(size_t readID);

    std::string filename;
};


#endif //BSG_LONGREADSDATASTORE_HPP
