//
// Created by Luis Yanes (EI) on 23/03/2018.
//

#include "LongReadsDatastore.hpp"

void LongReadsDatastore::load_index(std::string &file) {
    filename = file;
    std::ifstream ifs(file, std::ios_base::binary);
    std::getline(ifs, long_read_file, '\0');
    long_read_fd = ::open(long_read_file.c_str(), O_RDONLY);
    if (long_read_fd == -1) {
        perror("fstat");
        exit(EXIT_FAILURE);
    }
    read_rtfr(ifs);

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

uint32_t LongReadsDatastore::build_from_fastq(std::string long_read_file) {
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
    rc = posix_memalign((void **) &buffer, page_size, page_size);
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
        } while (ptr < sz_read);
    }
    free(buffer);
    return static_cast<uint32_t>(read_to_fileRecord.size());
}

void LongReadsDatastore::read_rtfr(std::ifstream &ifs) {
    uint64_t s;
    ifs.read((char *) &s, sizeof(s));
    read_to_fileRecord.resize(s);
    ifs.read((char *) read_to_fileRecord.data(), s*sizeof(read_to_fileRecord[0]));
}

void LongReadsDatastore::write_rtfr(std::ofstream &output_file) {
    //read filename
    uint64_t s=read_to_fileRecord.size();
    output_file.write((char *) &s,sizeof(s));
    output_file.write((char *)read_to_fileRecord.data(),s*sizeof(read_to_fileRecord[0]));
}

void LongReadsDatastore::read(std::ifstream &ifs) {
    uint64_t s;
    ifs.read((char *) &s, sizeof(s));
    filename.resize(s);
    ifs.read((char *) filename.data(), s*sizeof(filename[0]));
    load_index(filename);
}

void LongReadsDatastore::write(std::ofstream &output_file) {
    //read filename
    uint64_t s=filename.size();
    output_file.write((char *) &s,sizeof(s));
    output_file.write((char *)filename.data(),s*sizeof(filename[0]));
}

std::string LongReadsDatastore::get_read_sequence(size_t readID) {
    // Calculate the number of pages to mmap
    auto fileOffset = read_to_fileRecord[readID].offset+1;
    auto recordSize = read_to_fileRecord[readID].record_size;
    // Copy out the part where the sequence by
    std::string result((char*)mapped_readfile+fileOffset, recordSize);
    return result;
}

LongReadsDatastore::LongReadsDatastore(std::string filename) {
    load_index(filename);
}

LongReadsDatastore::LongReadsDatastore(std::string long_read_file, std::string output_file) {
    filename = output_file;
    build_from_fastq(long_read_file);
    std::ofstream ofs(output_file, std::ios_base::out | std::ios_base::trunc | std::ios_base::binary);
    ofs.write(long_read_file.data(), long_read_file.size()+1);
    write_rtfr(ofs);

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
