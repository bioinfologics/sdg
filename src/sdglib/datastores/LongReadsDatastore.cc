#include <memory>

//
// Created by Luis Yanes (EI) on 23/03/2018.
//

#include "LongReadsDatastore.hpp"
#include <sdglib/utilities/OutputLog.hpp>
#include <sdglib/mappers/LongReadsMapper.hpp>
#include <sdglib/workspace/WorkSpace.hpp>
#include <algorithm>

const sdgVersion_t LongReadsDatastore::min_compat = 0x0003;

LongReadsDatastore::LongReadsDatastore(WorkSpace &ws, std::string filename) : filename(filename), ws(ws), mapper(ws, *this) {
    load_index(filename);
}

LongReadsDatastore::LongReadsDatastore(WorkSpace &ws, const std::string &long_read_file, const std::string &output_file) : filename(output_file), ws(ws), mapper(ws, *this) {
    uint32_t nReads(1);
    std::ofstream ofs(output_file, std::ios_base::out | std::ios_base::trunc | std::ios_base::binary);
    if (!ofs) {
        std::cerr << "Failed to open " << output_file <<": " << strerror(errno);
        throw std::runtime_error("Could not open " + output_file);
    }
    std::streampos fPos;
    ofs.write((const char *) &SDG_MAGIC, sizeof(SDG_MAGIC));
    ofs.write((const char *) &SDG_VN, sizeof(SDG_VN));
    SDG_FILETYPE type(LongDS_FT);
    ofs.write((char *) &type, sizeof(type));

    ofs.write((char*) &nReads, sizeof(nReads));
    ofs.write((char*) &fPos, sizeof(fPos));
    nReads = dump_seqs_create_index(ofs, long_read_file); // Build read_to_fileRecord
    fPos = ofs.tellp();                             // Write position after reads
    sdglib::write_flat_vector(ofs, read_to_fileRecord);
    ofs.seekp(sizeof(SDG_MAGIC)+sizeof(SDG_VN)+sizeof(type));                                   // Go to top and dump # reads and position of index
    ofs.write((char*) &nReads, sizeof(nReads));     // Dump # of reads
    ofs.write((char*) &fPos, sizeof(fPos));         // Dump index
    ofs.flush();                                    // Make sure everything has been written
    sdglib::OutputLog(sdglib::LogLevels::INFO)<<"Built datastore with "<<size()<<" reads"<<std::endl;
    fd = open(filename.data(), O_RDONLY);
    if (!fd) {
        std::cerr << "Failed to open " << filename << ": " << strerror(errno);
        throw std::runtime_error("Could not open " + filename);
    }

}

LongReadsDatastore::LongReadsDatastore(WorkSpace &ws, std::ifstream &infile) : filename(), ws(ws), mapper(ws, *this) {
    read(infile);
    mapper.read(infile);
}

LongReadsDatastore::LongReadsDatastore(WorkSpace &ws, std::string default_name, const std::string &filename, std::ifstream &input_file) : filename(filename), ws(ws), mapper(ws, *this) {
    sdgMagic_t magic;
    sdgVersion_t version;
    SDG_FILETYPE type;
    input_file.read((char *) &magic, sizeof(magic));
    input_file.read((char *) &version, sizeof(version));
    input_file.read((char *) &type, sizeof(type));

    if (magic != SDG_MAGIC) {
        throw std::runtime_error(filename + " appears to be corrupted");
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
    std::streampos fPos;
    input_file.read((char*) &fPos, sizeof(fPos));

    sdglib::read_string(input_file, default_name);

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
    fd = open(filename.data(), O_RDONLY);
    if (!fd) {
        std::cerr << "Failed to open " << filename << ": " << strerror(errno);
        throw std::runtime_error("Could not open " + filename);
    }

}

LongReadsDatastore::LongReadsDatastore(WorkSpace &ws, LongReadsDatastore &o) : filename(o.filename), ws(ws), mapper(ws, *this) {
    read_to_fileRecord = o.read_to_fileRecord;
    name = o.name;
    default_name = o.default_name;
    mapper.first_mapping = o.mapper.first_mapping;
    mapper.read_paths = o.mapper.read_paths;
    mapper.all_paths_between = o.mapper.all_paths_between;
    fd = open(filename.data(), O_RDONLY);
    if (!fd) {
        std::cerr << "Failed to open " << filename << ": " << strerror(errno);
        throw std::runtime_error("Could not open " + filename);
    }
}

LongReadsDatastore &LongReadsDatastore::operator=(LongReadsDatastore const &o) {
    if (&o == this) return *this;

    filename = o.filename;
    ws = o.ws;
    read_to_fileRecord = o.read_to_fileRecord;

    mapper = o.mapper;
    default_name = o.default_name;
    name = o.name;
    fd = open(filename.data(), O_RDONLY);
    if (!fd) {
        std::cerr << "Failed to open " << filename <<": " << strerror(errno);
        throw std::runtime_error("Could not open " + filename + " when using the operator= of LongReadsDatastore");
    }
    return *this;
}

LongReadsDatastore::LongReadsDatastore(const LongReadsDatastore &o) :
        ws(o.ws),
        mapper(*this, o.mapper),
        read_to_fileRecord(o.read_to_fileRecord),
        name(o.name),
        default_name(o.default_name),
        filename(o.filename)
{
    fd = open(filename.data(), O_RDONLY);
}

void LongReadsDatastore::load_index(std::string &file) {
    filename = file;
    fd = open(filename.data(), O_RDONLY);
    std::ifstream input_file(file, std::ios_base::binary);
    if (!input_file) {
        std::cerr << "Failed to open " << file <<": " << strerror(errno);
        throw std::runtime_error("Could not open " + file);
    }

    uint64_t nReads(0);
    std::streampos fPos;

    sdgMagic_t magic;
    sdgVersion_t version;
    SDG_FILETYPE type;
    input_file.read((char *) &magic, sizeof(magic));
    input_file.read((char *) &version, sizeof(version));
    input_file.read((char *) &type, sizeof(type));

    if (magic != SDG_MAGIC) {
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

    sdglib::read_string(input_file, default_name);

    input_file.seekg(fPos);
    sdglib::read_flat_vector(input_file, read_to_fileRecord);

    sdglib::OutputLog()<<"LongReadsDatastore open: "<<filename<<" Total reads: " <<size()-1<<std::endl;
}

void LongReadsDatastore::build_from_fastq(const std::string &output_file, const std::string &default_name, const std::string &long_read_file, size_t min_size) {
    uint64_t nReads(1);
    std::vector<ReadPosSize> read_to_file_record{ReadPosSize{0,0}};
    std::ofstream ofs(output_file, std::ios_base::out | std::ios_base::trunc | std::ios_base::binary);
    if (!ofs) {
        std::cerr << "Failed to open output in " << output_file <<": " << strerror(errno);
        throw std::runtime_error("Could not open " + output_file);
    }
    std::streampos fPos;
    ofs.write((const char *) &SDG_MAGIC, sizeof(SDG_MAGIC));
    ofs.write((const char *) &SDG_VN, sizeof(min_compat));
    SDG_FILETYPE type(LongDS_FT);
    ofs.write((char *) &type, sizeof(type));

    ofs.write((char*) &nReads, sizeof(nReads));
    ofs.write((char*) &fPos, sizeof(fPos));

    sdglib::write_string(ofs, default_name);

    std::ifstream fastq_ifstream(long_read_file);
    if (!fastq_ifstream) {
        std::cerr << "Failed to open input in " << long_read_file << ": " << strerror(errno);
        throw std::runtime_error("Could not open " + long_read_file);
    }
    FastqReader<FastqRecord> reader({}, long_read_file);
    FastqRecord rec;
    while(reader.next_record(rec)) {
        if (!rec.seq.empty() and (min_size==0 or rec.seq.size()>=min_size)) {
            uint32_t size = rec.seq.size();
            auto offset = ofs.tellp();
            read_to_file_record.emplace_back((off_t)offset,size);
            ofs.write((char*)rec.seq.c_str(), size+1);//+1 writes the \0
        }
        ++nReads;
    }
    fPos = ofs.tellp();                             // Write position after reads
    sdglib::write_flat_vector(ofs, read_to_file_record);

    // Write empty mappings
    uint8_t k(0);
    ofs.write(reinterpret_cast<const char *>(&k), sizeof(k));
    std::vector<LongReadMapping> mappings;
    sdglib::write_flat_vector(ofs, mappings);

    ofs.seekp(sizeof(SDG_MAGIC)+sizeof(SDG_VN)+sizeof(type));                                   // Go to top and dump # reads and position of index
    ofs.write((char*) &nReads, sizeof(nReads));     // Dump # of reads
    ofs.write((char*) &fPos, sizeof(fPos));         // Dump index

    ofs.flush();                                    // Make sure everything has been written
    sdglib::OutputLog(sdglib::LogLevels::INFO)<<"Built datastore with "<<read_to_file_record.size()-1<<" reads"<<std::endl;
}

uint32_t LongReadsDatastore::dump_seqs_create_index(std::ofstream &outf, const std::string &long_read_file) {
    // open the file
    std::ifstream fastq_ifstream(long_read_file);
    if (!fastq_ifstream) {
        std::cerr << "Failed to open " << long_read_file <<": " << strerror(errno);
        throw std::runtime_error("Could not open " + long_read_file);
    }
    FastqReader<FastqRecord> reader({}, long_read_file);
    FastqRecord rec;
    while(reader.next_record(rec)) {
        if (!rec.seq.empty()) {
            uint32_t size = rec.seq.size();
            auto offset = outf.tellp();
            read_to_fileRecord.emplace_back((off_t)offset,size);
            outf.write((char*)rec.seq.c_str(), size+1);//+1 writes the \0
        }
    }
    return static_cast<uint32_t>(read_to_fileRecord.size());
}

void LongReadsDatastore::print_status() const {
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
    sdglib::read_string(ifs, filename);
    sdglib::read_string(ifs, name);

    load_index(filename);
    fd = open(filename.data(), O_RDONLY);
    if (!fd) {
        std::cerr << "Failed to open " << filename << ": " << strerror(errno);
        throw std::runtime_error("Could not open " + filename);
    }
}

void LongReadsDatastore::write(std::ofstream &output_file) {
    //read filename
    sdglib::write_string(output_file, filename);
    sdglib::write_string(output_file, name);

    mapper.write(output_file);
}

std::string LongReadsDatastore::get_read_sequence(size_t readID) const {
    char buffer[read_to_fileRecord[readID].record_size+1];
    size_t read_offset_in_file=read_to_fileRecord[readID].offset;
    ::lseek(fd,read_offset_in_file,SEEK_SET);
    ::read(fd, buffer,sizeof(buffer));
    return std::string(buffer);
}

void LongReadsDatastore::write_selection(std::ofstream &output_file, const std::vector<uint64_t> &read_ids) {
    unsigned long size(read_ids.size());
    auto rsb=ReadSequenceBuffer(*this);
    output_file.write((const char *) &SDG_MAGIC, sizeof(SDG_MAGIC));
    output_file.write((const char *) &SDG_VN, sizeof(SDG_VN));
    SDG_FILETYPE type(LongDS_FT);
    output_file.write((char *) &type, sizeof(type));

    output_file.write((char *) &size, sizeof(size)); // How many reads we will write in the file
    const char* seq_ptr;
    for (unsigned long long read_id : read_ids) {
        seq_ptr = rsb.get_read_sequence(read_id);
        std::string seq(seq_ptr);
        size=seq.size();
        output_file.write((char *) &size,sizeof(size)); // Read length
        output_file.write(seq.data(),seq.size());       // Read sequence
    }
}

LongReadsDatastore::~LongReadsDatastore() {
    close(fd);
}

std::ostream &operator<<(std::ostream &os, const LongReadsDatastore &lords) {
    os << "LongReadsDatastore" << std::endl;
    os << "Name: " << ( (lords.name.empty()) ? lords.default_name : lords.name) << std::endl;
    lords.print_status();
}

