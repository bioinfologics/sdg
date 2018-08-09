//
// Created by Luis Yanes (EI) on 23/03/2018.
//

#include <sglib/logger/OutputLog.h>
#include "LongReadsDatastore.hpp"

void LongReadsDatastore::load_index(std::string &file) {
    filename = file;
    std::ifstream ifs(file, std::ios_base::binary);
    std::getline(ifs, file_containing_long_read_sequence, '\0');
    lr_sequence_fd = ::open(file_containing_long_read_sequence.c_str(), O_RDONLY);
    if (lr_sequence_fd == -1) {
        perror("fstat");
        exit(EXIT_FAILURE);
    }
    read_rtfr(ifs);

    struct stat sb;

    if (fstat(lr_sequence_fd, &sb) == -1) {         /* To obtain file size */
        perror("fstat");
        exit(EXIT_FAILURE);
    }

    auto length = sb.st_size - offset;
    off_t pa_offset = offset & (~(page_size-1));

    this->mapped_readfile = ::mmap(NULL, length + offset - pa_offset, PROT_READ,
                                   MAP_PRIVATE, lr_sequence_fd, pa_offset);
    sglib::OutputLog()<<"LongReadsDatastore open: "<<filename<<" Total reads: " <<size()<<std::endl;
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

std::string LongReadsDatastore::get_read_sequence(size_t readID) const {
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

    lr_sequence_fd = ::open(long_read_file.c_str(), O_RDONLY);
    if (lr_sequence_fd == -1) {
        perror("fstat");
        exit(EXIT_FAILURE);
    }
    struct stat sb;

    if (fstat(lr_sequence_fd, &sb) == -1) {         /* To obtain file size */
        perror("fstat");
        exit(EXIT_FAILURE);
    }

    auto length = sb.st_size - offset;
    off_t pa_offset = offset & (~(page_size-1));

    this->mapped_readfile = ::mmap(NULL, length + offset - pa_offset, PROT_READ,
                                   MAP_PRIVATE, lr_sequence_fd, pa_offset);

}

void LongReadsDatastore::write_selection(std::ofstream &output_file, const std::vector<uint64_t> &read_ids) const {
    uint32_t size(this->file_containing_long_read_sequence.size());


    size = read_ids.size();
    output_file.write((char *) &size, sizeof(size)); // How many reads we will write in the file
    for (auto i=0;i<read_ids.size();i++) {
        std::string seq(get_read_sequence(read_ids[i]));
        size=seq.size();
        output_file.write((char *) &size,sizeof(size)); // Read length
        output_file.write(seq.data(),seq.size());       // Read sequence
    }
}

void LongReadsDatastore::load_from_stream(std::string file_name, std::ifstream &input_file) {
    file_containing_long_read_sequence = file_name;

    lr_sequence_fd = ::open(file_containing_long_read_sequence.c_str(), O_RDONLY);
    if (lr_sequence_fd == -1) {
        perror("fstat");
        exit(EXIT_FAILURE);
    }

    // Equivalente de read_rtfr
    size_t numReads(0);
    input_file.read((char *) numReads, sizeof(numReads));
    read_to_fileRecord.resize(numReads);
    for (auto i=0; i<numReads; i++) {
        read_to_fileRecord[i].offset=input_file.tellg();
        uint32_t readSize(0);
        input_file.read((char *) &readSize, sizeof(readSize));
        std::string seq;
        seq.resize(readSize);
        input_file.read(&seq[0], readSize);
        read_to_fileRecord[i].record_size=readSize;
    }

    struct stat sb;

    if (fstat(lr_sequence_fd, &sb) == -1) {         /* To obtain file size */
        perror("fstat");
        exit(EXIT_FAILURE);
    }

    auto length = sb.st_size - offset;
    off_t pa_offset = offset & (~(page_size-1));

    this->mapped_readfile = ::mmap(NULL, length + offset - pa_offset, PROT_READ,
                                   MAP_PRIVATE, lr_sequence_fd, pa_offset);

}
