//
// Created by Luis Yanes (EI) on 23/03/2018.
//

#include "LongReadsDatastore.hpp"
#include <sdglib/logger/OutputLog.hpp>
#include <sdglib/mappers/LongReadsMapper.hpp>
#include <sdglib/workspace/WorkSpace.hpp>
#include <algorithm>

const bsgVersion_t LongReadsDatastore::min_compat = 0x0002;

void LongReadsDatastore::load_index(std::string &file) {
    filename = file;

    std::ifstream input_file(file, std::ios_base::binary);
    if (!input_file) {
        std::cerr << "Failed to open " << file <<": " << strerror(errno);
        throw std::runtime_error("Could not open " + file);

    }

    uint64_t nReads(0);
    std::streampos fPos;

    bsgMagic_t magic;
    bsgVersion_t version;
    BSG_FILETYPE type;
    input_file.read((char *) &magic, sizeof(magic));
    input_file.read((char *) &version, sizeof(version));
    input_file.read((char *) &type, sizeof(type));

    if (magic != BSG_MAGIC) {
        throw std::runtime_error("This file appears to be corrupted: " + file);
    }

    if (version < min_compat) {
        throw std::runtime_error("Incompatible version" + std::to_string(version));
    }

    if (type != LongDS_FT) {
        throw std::runtime_error("Incompatible file type" + std::to_string(type));
    }

    input_file.read((char*)&nReads, sizeof(nReads));

    input_file.read((char*)&fPos, sizeof(fPos));
    input_file.seekg(fPos);
    read_to_fileRecord.resize(nReads);
    input_file.read((char *) read_to_fileRecord.data(), read_to_fileRecord.size()*sizeof(read_to_fileRecord[0]));

    sdglib::OutputLog()<<"LongReadsDatastore open: "<<filename<<" Total reads: " <<size()-1<<std::endl;
}

void LongReadsDatastore::build_from_fastq(std::string output_file, std::string long_read_file) {
    uint64_t nReads(0);
    std::vector<ReadPosSize> read_to_file_record;
    std::ofstream ofs(output_file, std::ios_base::out | std::ios_base::trunc | std::ios_base::binary);
    if (!ofs) {
        std::cerr << "Failed to open output in " << output_file <<": " << strerror(errno);
        throw std::runtime_error("Could not open " + output_file);
    }
    std::streampos fPos;
    ofs.write((const char *) &BSG_MAGIC, sizeof(BSG_MAGIC));
    ofs.write((const char *) &BSG_VN, sizeof(BSG_VN));
    BSG_FILETYPE type(LongDS_FT);
    ofs.write((char *) &type, sizeof(type));

    ofs.write((char*) &nReads, sizeof(nReads));
    ofs.write((char*) &fPos, sizeof(fPos));
    std::ifstream fastq_ifstream(long_read_file);
    if (!fastq_ifstream) {
        std::cerr << "Failed to open input in " << long_read_file << ": " << strerror(errno);
        throw std::runtime_error("Could not open " + long_read_file);
    }
    std::string name,seq,p,qual;
    while(fastq_ifstream.good()) {
        std::getline(fastq_ifstream, name);
        std::getline(fastq_ifstream, seq);
        if (!seq.empty()) {
            uint32_t size = seq.size();
            auto offset = ofs.tellp();
            read_to_file_record.emplace_back((off_t)offset,size);
            ofs.write((char*)seq.c_str(), size+1);//+1 writes the \0
        }
        std::getline(fastq_ifstream, p);
        std::getline(fastq_ifstream, qual);
        ++nReads;
    }
    read_to_file_record.pop_back();
    fPos = ofs.tellp();                             // Write position after reads
    ofs.write((char *)read_to_file_record.data(),read_to_file_record.size()*sizeof(read_to_fileRecord[0]));
    ofs.seekp(sizeof(BSG_MAGIC)+sizeof(BSG_VN)+sizeof(type));                                   // Go to top and dump # reads and position of index
    ofs.write((char*) &nReads, sizeof(nReads));     // Dump # of reads
    ofs.write((char*) &fPos, sizeof(fPos));         // Dump index
    ofs.flush();                                    // Make sure everything has been written
    sdglib::OutputLog(sdglib::LogLevels::INFO)<<"Built datastore with "<<read_to_file_record.size()-1<<" reads"<<std::endl;
}
uint32_t LongReadsDatastore::build_from_fastq(std::ofstream &outf, std::string long_read_file) {
    // open the file
    uint32_t numRecords(0);
    std::ifstream fastq_ifstream(long_read_file);
    if (!fastq_ifstream) {
        std::cerr << "Failed to open " << long_read_file <<": " << strerror(errno);
        throw std::runtime_error("Could not open " + long_read_file);
    }
    std::string name,seq,p,qual;
    while(fastq_ifstream.good()) {
        std::getline(fastq_ifstream, name);
        std::getline(fastq_ifstream, seq);
        if (!seq.empty()) {
            uint32_t size = seq.size();
            auto offset = outf.tellp();
            read_to_fileRecord.emplace_back((off_t)offset,size);
            outf.write((char*)seq.c_str(), size+1);//+1 writes the \0
        }
        std::getline(fastq_ifstream, p);
        std::getline(fastq_ifstream, qual);
        ++numRecords;
    }
    read_to_fileRecord.pop_back();
    return static_cast<uint32_t>(read_to_fileRecord.size());
}

void LongReadsDatastore::print_status() {
    auto &log_no_date=sdglib::OutputLog(sdglib::LogLevels::INFO,false);
    std::vector<uint64_t> read_sizes;
    uint64_t total_size=0;
    for (auto ri=1;ri<size();++ri) {
        auto s=read_to_fileRecord[ri].record_size;
            total_size += s;
            read_sizes.push_back(s);
    }
    std::sort(read_sizes.rbegin(),read_sizes.rend());
    auto on20s=total_size*.2;
    auto on50s=total_size*.5;
    auto on80s=total_size*.8;
    sdglib::OutputLog() <<"The datastore on " << filename << " contains "<<read_sizes.size()<<" reads with "<<total_size<<"bp, ";
    uint64_t acc=0;
    for (auto s:read_sizes){
        if (acc==0)  log_no_date<<"N0: "<<s<<"bp  ";
        if (acc<on20s and acc+s>on20s) log_no_date<<"N20: "<<s<<"bp  ";
        if (acc<on50s and acc+s>on50s) log_no_date<<"N50: "<<s<<"bp  ";
        if (acc<on80s and acc+s>on80s) log_no_date<<"N80: "<<s<<"bp  ";
        acc+=s;
        if (acc==total_size)  log_no_date<<"N100: "<<s<<"bp  ";
    }
    log_no_date<<std::endl;

    mapper.print_status();
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

LongReadsDatastore::LongReadsDatastore(WorkSpace &ws, std::string filename) : ws(ws), mapper(ws, *this) {
    load_index(filename);
}

LongReadsDatastore::LongReadsDatastore(WorkSpace &ws, std::string long_read_file, std::string output_file) : ws(ws), mapper(ws, *this) {
    filename = output_file;

    uint32_t nReads(0);
    std::ofstream ofs(output_file, std::ios_base::out | std::ios_base::trunc | std::ios_base::binary);
    if (!ofs) {
        std::cerr << "Failed to open " << output_file <<": " << strerror(errno);
        throw std::runtime_error("Could not open " + output_file);
    }
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
    sdglib::OutputLog(sdglib::LogLevels::INFO)<<"Built datastore with "<<size()-1<<" reads"<<std::endl;
}

LongReadsDatastore::LongReadsDatastore(WorkSpace &ws, std::ifstream &infile) : ws(ws), mapper(ws, *this) {
    read(infile);
}

LongReadsDatastore::LongReadsDatastore(WorkSpace &ws, std::string file_name, std::ifstream &input_file) : ws(ws), mapper(ws, *this) {
    file_containing_long_read_sequence = file_name;

    bsgMagic_t magic;
    bsgVersion_t version;
    BSG_FILETYPE type;
    input_file.read((char *) &magic, sizeof(magic));
    input_file.read((char *) &version, sizeof(version));
    input_file.read((char *) &type, sizeof(type));

    if (magic != BSG_MAGIC) {
        throw std::runtime_error(file_name + " appears to be corrupted");
    }

    if (version < min_compat) {
        throw std::runtime_error("Incompatible version");
    }

    if (type != LongDS_FT) {
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

LongReadsDatastore::LongReadsDatastore(WorkSpace &ws, LongReadsDatastore &o) : ws(ws), mapper(ws, *this) {
    this->filename = o.filename;
    this->file_containing_long_read_sequence = o.file_containing_long_read_sequence;
    this->read_to_fileRecord = o.read_to_fileRecord;

    this->mapper.reads_in_node = o.mapper.reads_in_node;
    this->mapper.filtered_read_mappings = o.mapper.filtered_read_mappings;
    this->mapper.first_mapping = o.mapper.first_mapping;
    this->mapper.read_paths = o.mapper.read_paths;
    this->mapper.all_paths_between = o.mapper.all_paths_between;
}

LongReadsDatastore &LongReadsDatastore::operator=(LongReadsDatastore const &o) {
    if (&o == this) return *this;

    this->filename = o.filename;
    this->ws = o.ws;
    this->read_to_fileRecord = o.read_to_fileRecord;
    this->file_containing_long_read_sequence = o.file_containing_long_read_sequence;

    this->mapper = o.mapper;
    return *this;
}

LongReadsDatastore::LongReadsDatastore(const LongReadsDatastore &o) :
    ws(o.ws),
    mapper(*this, o.mapper),
    read_to_fileRecord(o.read_to_fileRecord),
    filename(o.filename),
    file_containing_long_read_sequence(o.file_containing_long_read_sequence)
    {}

void BufferedSequenceGetter::write_selection(std::ofstream &output_file, const std::vector<uint64_t> &read_ids) {
    unsigned long size(read_ids.size());

    output_file.write((const char *) &BSG_MAGIC, sizeof(BSG_MAGIC));
    output_file.write((const char *) &BSG_VN, sizeof(BSG_VN));
    BSG_FILETYPE type(LongDS_FT);
    output_file.write((char *) &type, sizeof(type));

    output_file.write((char *) &size, sizeof(size)); // How many reads we will write in the file
    const char* seq_ptr;
    for (auto i=0;i<read_ids.size();i++) {
        seq_ptr = get_read_sequence(read_ids[i]);
        std::string seq(seq_ptr);
        size=seq.size();
        output_file.write((char *) &size,sizeof(size)); // Read length
        output_file.write(seq.data(),seq.size());       // Read sequence
    }
}

const char * BufferedSequenceGetter::get_read_sequence(uint64_t readID) {
    off_t read_offset_in_file=datastore.read_to_fileRecord[readID].offset;
    if (chunk_size < datastore.read_to_fileRecord[readID].record_size) {
        throw std::runtime_error(
                "Reading from " + this->datastore.filename +
                " failed!\nThe size of the buffer chunk is smaller than read " +
                std::to_string(readID) + " increase the chunk_size so this read fits");
    }
    if (read_offset_in_file<buffer_offset or read_offset_in_file+chunk_size>buffer_offset+bufsize) {
        buffer_offset=read_offset_in_file;
        lseek(fd,read_offset_in_file,SEEK_SET);
        size_t sz_read=0, total_left = 0;
        total_left = std::min((uint64_t) bufsize, (uint64_t) total_size-read_offset_in_file);
        char * buffer_pointer = buffer;
        while(total_left>0) {
            ssize_t current = read(fd, buffer_pointer, total_left);
            if (current < 0) {
                throw std::runtime_error("Reading from " + this->datastore.filename + " failed!");
            } else {
                sz_read += current;
                total_left -= current;
                buffer_pointer += current;
            }
        }
    }
    return buffer+(read_offset_in_file-buffer_offset);
}
