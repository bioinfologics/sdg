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

public:
    std::vector< ReadPosSize > read_to_fileRecord{ReadPosSize(0,0)};
    LongReadsDatastore() = delete;

    /**
     * Initialize from already created index
     * @param filename
     */
    LongReadsDatastore(std::string filename) {
        load_index(filename);
    }

    /**
     * Initialize from long_read_file then store the index
     * @param long_read_file
     * @param output_file
     */
    LongReadsDatastore(std::string long_read_file, std::string output_file) {
        build_from_fastq(long_read_file, output_file);
        std::ofstream ofs(output_file, std::ios_base::out | std::ios_base::trunc | std::ios_base::binary);
        ofs.write(long_read_file.data(), long_read_file.size()+1);
        write(ofs);

        long_read_fd = ::open(long_read_file.c_str(), O_RDONLY);
        if (long_read_fd == -1) {
            perror("fstat");
            exit(EXIT_FAILURE);
        }
        struct stat sb;

        if (fstat(long_read_fd, &sb) == -1) {         /* To obtain file size */
            perror("fstat");
            exit(EXIT_FAILURE);
        }

        auto length = sb.st_size - offset;
        off_t pa_offset = offset & (~(page_size-1));

        this->mapped_readfile = ::mmap(NULL, length + offset - pa_offset, PROT_READ,
                                       MAP_PRIVATE, long_read_fd, pa_offset);

    }

    void load_index(std::string &file) {
        std::ifstream ifs(file, std::ios_base::binary);
        std::getline(ifs, long_read_file, '\0');
        long_read_fd = ::open(long_read_file.c_str(), O_RDONLY);
        if (long_read_fd == -1) {
            perror("fstat");
            exit(EXIT_FAILURE);
        }
        read(ifs);

        struct stat sb;

        if (fstat(long_read_fd, &sb) == -1) {         /* To obtain file size */
            perror("fstat");
            exit(EXIT_FAILURE);
        }

        auto length = sb.st_size - offset;
        off_t pa_offset = offset & (~(page_size-1));

        this->mapped_readfile = ::mmap(NULL, length + offset - pa_offset, PROT_READ,
                    MAP_PRIVATE, long_read_fd, pa_offset);

    }

    uint32_t build_from_fastq(std::string long_read_file, std::string output_file) {
        // open the file
        read_to_fileRecord.reserve(100000);
        uint32_t numRecords(0);
        int infile = ::open(long_read_file.c_str(), O_RDONLY);
        if (!infile) {
            return 0;
        }
        // read each record
        char *buffer;
        int rc;
        rc = posix_memalign((void **) &buffer, page_size, 1);
        if (0 != rc) {
            return 0;
        }

        //res.reserve(100000);
        ssize_t sz=0, sz_read=1;
        uint32_t eolCounter = 0;
        uint32_t readSize = 0;
        off_t filePos = 0;
        while (sz_read > 0
               && (sz_read = ::read(infile, buffer, page_size)) > 0) {
            size_t ptr = 0;
            do {
                if (buffer[ptr] == '\n'){
                    eolCounter++;
                    auto modEol(eolCounter % 4);
                    if (modEol == 1) {
                        read_to_fileRecord.push_back(ReadPosSize(filePos,0));
                        readSize=0;
                    } else if (modEol == 2) {
                        read_to_fileRecord.back().record_size = readSize;
                    }
                }
                ++readSize;
                ++filePos;
                ++ptr;
            } while (ptr != sz_read);
        }
        free(buffer);
        return static_cast<uint32_t>(read_to_fileRecord.size());
    }

    void read(std::ifstream &ifs) {
        uint64_t s;
        ifs.read((char *) &s, sizeof(s));
        read_to_fileRecord.resize(s);
        ifs.read((char *) read_to_fileRecord.data(), s*sizeof(read_to_fileRecord[0]));
    }
    void write(std::ofstream &output_file){
        uint64_t s=read_to_fileRecord.size();
        output_file.write((char *) &s,sizeof(s));
        output_file.write((char *)read_to_fileRecord.data(),s*sizeof(read_to_fileRecord[0]));
    }

    // MMap a section of the file then get the mmapped pages
    std::string get_read_sequence(size_t readID){
        // Calculate the number of pages to mmap


        // Copy out the part where the sequence by
        std::string result((char*)mapped_readfile+read_to_fileRecord[readID].offset+1, read_to_fileRecord[readID].record_size-1);
        ::munmap(this->mapped_readfile, read_to_fileRecord[readID].record_size);
        return result;
    }
};


#endif //BSG_LONGREADSDATASTORE_HPP
