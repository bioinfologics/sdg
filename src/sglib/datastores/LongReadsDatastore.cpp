//
// Created by Luis Yanes (EI) on 23/03/2018.
//

#include <sglib/logger/OutputLog.h>
#include "LongReadsDatastore.hpp"

const bsgVersion_t LongReadsDatastore::min_compat = 0x0001;

void LongReadsDatastore::load_index(std::string &file) {
    filename = file;

    std::ifstream input_file(file, std::ios_base::binary);

    uint32_t nReads(0);
    std::streampos fPos;

    bsgMagic_t magic;
    bsgVersion_t version;
    BSG_FILETYPE type;
    input_file.read((char *) &magic, sizeof(magic));
    input_file.read((char *) &version, sizeof(version));
    input_file.read((char *) &type, sizeof(type));

    if (magic != BSG_MAGIC) {
        std::cerr << "This file seems to be corrupted" << std::endl;
        throw std::runtime_error("This file appears to be corrupted");
    }

    if (version < min_compat) {
        std::cerr << "This version of the file is not compatible with the current build, please update" << std::endl;
        throw std::runtime_error("Incompatible version");
    }

    if (type != LongDS_FT) {
        std::cerr << "This file is not compatible with this type" << std::endl;
        throw std::runtime_error("Incompatible file type");
    }

    input_file.read((char*)&nReads, sizeof(nReads));
    input_file.read((char*)&fPos, sizeof(fPos));
    input_file.seekg(fPos);
    read_to_fileRecord.resize(nReads);
    input_file.read((char *) read_to_fileRecord.data(), read_to_fileRecord.size()*sizeof(read_to_fileRecord[0]));

    sglib::OutputLog()<<"LongReadsDatastore open: "<<filename<<" Total reads: " <<size()-1<<std::endl;
}

uint32_t LongReadsDatastore::build_from_fastq(std::ofstream &outf, std::string long_read_file) {
    // open the file
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
            auto offset = outf.tellp();
            outf.write((char*)&size, sizeof(size));
            read_to_fileRecord.emplace_back((off_t)offset+sizeof(BSG_MAGIC)+ sizeof(BSG_VN)+sizeof(BSG_FILETYPE),size);
            outf.write((char*)seq.c_str(), size);
        }
        std::getline(fastq_ifstream, p);
        std::getline(fastq_ifstream, qual);
        ++numRecords;
    }
    read_to_fileRecord.pop_back();
    return static_cast<uint32_t>(read_to_fileRecord.size());
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

LongReadsDatastore::LongReadsDatastore(std::string filename) {
    load_index(filename);
}

LongReadsDatastore::LongReadsDatastore(std::string long_read_file, std::string output_file) {
    filename = output_file;

    uint32_t nReads(0);
    std::ofstream ofs(output_file, std::ios_base::out | std::ios_base::trunc | std::ios_base::binary);
    std::streampos fPos;
    ofs.write((const char *) &BSG_MAGIC, sizeof(BSG_MAGIC));
    ofs.write((const char *) &BSG_VN, sizeof(BSG_VN));
    BSG_FILETYPE type(LongDS_FT);
    ofs.write((char *) &type, sizeof(type));

    ofs.write((char*) &nReads, sizeof(nReads));
    ofs.write((char*) &fPos, sizeof(fPos));
    nReads = build_from_fastq(ofs, long_read_file); // Build read_to_fileRecord
    fPos = ofs.tellp();                             // Write position after reads
    ofs.write((char *)read_to_fileRecord.data(),read_to_fileRecord.size()*sizeof(read_to_fileRecord[0]));
    ofs.seekp(sizeof(BSG_MAGIC)+sizeof(BSG_VN)+sizeof(type));                                   // Go to top and dump # reads and position of index
    ofs.write((char*) &nReads, sizeof(nReads));     // Dump # of reads
    ofs.write((char*) &fPos, sizeof(fPos));         // Dump index
    ofs.flush();                                    // Make sure everything has been written
    sglib::OutputLog(sglib::LogLevels::INFO)<<"Built datastore with "<<size()-1<<" reads"<<std::endl;
    lr_sequence_fd = ::open(output_file.c_str(), O_RDONLY);
    if (lr_sequence_fd == -1) {
        perror("fstat");
        throw std::runtime_error("Cannot open "+output_file);
    }

}

void LongReadsDatastore::load_from_stream(std::string file_name, std::ifstream &input_file) {
    file_containing_long_read_sequence = file_name;

    lr_sequence_fd = ::open(file_containing_long_read_sequence.c_str(), O_RDONLY);
    if (lr_sequence_fd == -1) {
        perror("fstat");
        throw std::runtime_error("Cannot open " + file_containing_long_read_sequence);
    }

    bsgMagic_t magic;
    bsgVersion_t version;
    BSG_FILETYPE type;
    input_file.read((char *) &magic, sizeof(magic));
    input_file.read((char *) &version, sizeof(version));
    input_file.read((char *) &type, sizeof(type));

    if (magic != BSG_MAGIC) {
        std::cerr << "This file seems to be corrupted" << std::endl;
        throw std::runtime_error("This file appears to be corrupted");
    }

    if (version < min_compat) {
        std::cerr << "This version of the file is not compatible with the current build, please update" << std::endl;
        throw std::runtime_error("Incompatible version");
    }

    if (type != LongDS_FT) {
        std::cerr << "This file is not compatible with this type" << std::endl;
        throw std::runtime_error("Incompatible file type");
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
}

void BufferedSequenceGetter::write_selection(std::ofstream &output_file, const std::vector<uint64_t> &read_ids) {
    unsigned long size(read_ids.size());

    output_file.write((const char *) &BSG_MAGIC, sizeof(BSG_MAGIC));
    output_file.write((const char *) &BSG_VN, sizeof(BSG_VN));
    BSG_FILETYPE type(LongDS_FT);
    output_file.write((char *) &type, sizeof(type));

    output_file.write((char *) &size, sizeof(size)); // How many reads we will write in the file
    for (auto i=0;i<read_ids.size();i++) {
        std::string seq(get_read_sequence(read_ids[i]));
        size=seq.size();
        output_file.write((char *) &size,sizeof(size)); // Read length
        output_file.write(seq.data(),seq.size());       // Read sequence
    }
}

std::string BufferedSequenceGetter::get_read_sequence(uint64_t readID) {
    off_t read_offset_in_file=datastore.read_to_fileRecord[readID].offset;
    if (read_offset_in_file<buffer_offset or read_offset_in_file+datastore.read_to_fileRecord[readID].record_size>buffer_offset+bufsize) {
        buffer_offset=read_offset_in_file;
        lseek(fd,read_offset_in_file,SEEK_SET);
        read(fd,buffer,bufsize);
    }
    return std::string(buffer+(read_offset_in_file-buffer_offset),buffer+(read_offset_in_file-buffer_offset)+datastore.read_to_fileRecord[readID].record_size);
}
