//
// Created by Luis Yanes (EI) on 23/03/2018.
//

#include <sglib/logger/OutputLog.h>
#include "LongReadsDatastore.hpp"

void LongReadsDatastore::load_index(std::string &file) {
    filename = file;

    std::ifstream ifs(file, std::ios_base::binary);

    uint64_t nReads(0);
    std::streampos fPos;

    ifs.read((char*)&nReads, sizeof(nReads));
    ifs.read((char*)&fPos, sizeof(fPos));
    ifs.seekg(fPos);
    read_rtfr(ifs);

    lr_sequence_fd = ::open(filename.c_str(), O_RDONLY);
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
    sglib::OutputLog()<<"LongReadsDatastore open: "<<filename<<" Total reads: " <<size()<<std::endl;
}

uint32_t LongReadsDatastore::build_from_fastq(std::ofstream &outf, std::string long_read_file) {
    // open the file
    read_to_fileRecord.reserve(100000);
    uint32_t numRecords(0);
    int infile = ::open(long_read_file.c_str(), O_RDONLY);
    if (!infile) {
        return 0;
    }
    std::ifstream fastq_ifstream(long_read_file);
    std::string name,seq,p,qual;
    while(fastq_ifstream.good()) {
        std::getline(fastq_ifstream, name);
        std::getline(fastq_ifstream, seq);
        if (!seq.empty()) {
            uint32_t size = seq.size();
            outf.write((char*)&size, sizeof(size));
            read_to_fileRecord.push_back(ReadPosSize(outf.tellp(),size));
            outf.write((char*)seq.c_str(), size);
        }
        std::getline(fastq_ifstream, p);
        std::getline(fastq_ifstream, qual);
        ++numRecords;
    }
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

    uint64_t nReads(0);
    std::ofstream ofs(output_file, std::ios_base::out | std::ios_base::trunc | std::ios_base::binary);
    std::streampos fPos;
    ofs.write((char*) &nReads, sizeof(nReads));
    ofs.write((char*) &fPos, sizeof(fPos));
    nReads = build_from_fastq(ofs, long_read_file); // Build read_to_fileRecord
    fPos = ofs.tellp();                             // Write position after reads
    write_rtfr(ofs);                                // Dump rtfr
    ofs.seekp(0);                                   // Go to top and dump # reads and position of index
    ofs.write((char*) &nReads, sizeof(nReads));     // Dump # of reads
    ofs.write((char*) &fPos, sizeof(fPos));         // Dump index
    ofs.flush();                                    // Make sure everything has been written
    lr_sequence_fd = ::open(output_file.c_str(), O_RDONLY);
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
