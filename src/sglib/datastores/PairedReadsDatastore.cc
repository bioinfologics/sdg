//
// Created by Bernardo Clavijo (EI) on 10/02/2018.
//

#include <sglib/logger/OutputLog.h>
#include <fstream>
#include <strings.h>
#include <cstring>
#include "PairedReadsDatastore.hpp"

void PairedReadsDatastore::build_from_fastq(std::string read1_filename,std::string read2_filename, std::string output_filename, int _min_rs, int _rs, size_t chunksize) {

    //std::cout<<"Memory used by every read's entry:"<< sizeof(PairedRead)<<std::endl;
    //read each read, put it on the index and on the appropriate tag
    readsize=_rs;
    sglib::OutputLog(sglib::LogLevels::INFO)<<"Creating Datastore from "<<read1_filename<<" | "<<read2_filename<<std::endl;
    auto fd1=fopen(read1_filename.c_str(),"r");
    auto fd2=fopen(read2_filename.c_str(),"r");
    char readbuffer[1000];
    //first, build an index of tags and offsets
    sglib::OutputLog()<<"Reading chunks of "<<chunksize<<" pairs"<<std::endl;
    std::vector<PairedReadData> readdatav;
    readdatav.reserve(chunksize);
    std::vector<std::ifstream> chunkfiles;
    std::vector<PairedReadData> next_in_chunk;
    PairedReadData currrent_read;
    //First, create the chunk files
    uint64_t pairs=0,discarded=0,truncated=0;
    std::ofstream output(output_filename.c_str());

    output.write((const char *) &readsize,sizeof(readsize));
    auto size_pos=output.tellp();
    output.write((const char *) &_size, sizeof(_size));//just to save the space!
    while (!feof(fd1) and !feof(fd2)) {

        if (NULL == fgets(readbuffer, 999, fd1)) continue;
        if (NULL == fgets(readbuffer, 999, fd1)) continue;
        currrent_read.seq1=std::string(readbuffer);
        if (NULL == fgets(readbuffer, 999, fd1)) continue;
        if (NULL == fgets(readbuffer, 999, fd1)) continue;
        if (currrent_read.seq1.back()=='\n') currrent_read.seq1.resize(currrent_read.seq1.size()-1);
        if (NULL == fgets(readbuffer, 999, fd2)) continue;
        if (NULL == fgets(readbuffer, 999, fd2)) continue;
        currrent_read.seq2=std::string(readbuffer);
        if (NULL == fgets(readbuffer, 999, fd2)) continue;
        if (NULL == fgets(readbuffer, 999, fd2)) continue;
        if (currrent_read.seq2.back()=='\n') currrent_read.seq2.resize(currrent_read.seq2.size()-1);
        if (currrent_read.seq1.size()<_min_rs or currrent_read.seq2.size()<_min_rs) {
            ++discarded;
            continue;
        }
        if (currrent_read.seq1.size()>_rs) {
            ++truncated;
            currrent_read.seq1.resize(_rs);
        }
        if (currrent_read.seq2.size()>_rs) {
            ++truncated;
            currrent_read.seq2.resize(_rs);
        }
        ++pairs;
        readdatav.push_back(currrent_read);
        if (readdatav.size()==chunksize){
            //dump
            char buffer[2*readsize+2];
            bzero(buffer,2*readsize+2);
            for (auto &r:readdatav){
                bzero(buffer,2*readsize+2);
                memcpy(buffer,r.seq1.data(),(r.seq1.size()>readsize ? readsize : r.seq1.size()));
                memcpy(buffer+readsize+1,r.seq2.data(),(r.seq2.size()>readsize ? readsize : r.seq2.size()));
                output.write(buffer,2*readsize+2);
            }
            sglib::OutputLog()<<readdatav.size()<<" pairs dumped..."<<std::endl;
            readdatav.clear();
        }
    }
    if (readdatav.size()>0) {
        //dump
        char buffer[2*readsize+2];
        bzero(buffer,2*readsize+2);
        for (auto &r:readdatav){
            bzero(buffer,2*readsize+2);
            memcpy(buffer,r.seq1.data(),(r.seq1.size()>readsize ? readsize : r.seq1.size()));
            memcpy(buffer+readsize+1,r.seq2.data(),(r.seq2.size()>readsize ? readsize : r.seq2.size()));
            output.write(buffer,2*readsize+2);
        }
        sglib::OutputLog()<<readdatav.size()<<" pairs dumped..."<<std::endl;
        readdatav.clear();
    }
    _size=pairs*2;
    output.seekp(size_pos);
    output.write((const char *) &_size, sizeof(_size));
    output.close();
    //DONE!
    sglib::OutputLog(sglib::LogLevels::INFO)<<discarded<<" pairs discarded due to short reads"<<std::endl;
    sglib::OutputLog(sglib::LogLevels::INFO)<<truncated<<" reads where truncated to "<<_rs<<"bp"<<std::endl;
    sglib::OutputLog(sglib::LogLevels::INFO)<<"Datastore with "<<_size<<" reads ("<<pairs<<" pairs)"<<std::endl;
    filename=output_filename;
    fd=fopen(filename.c_str(),"r");
}

void PairedReadsDatastore::read(std::ifstream &input_file) {
    //read filename
    uint64_t s;
    input_file.read((char *) &s, sizeof(s));
    filename.resize(s);
    input_file.read((char *) filename.data(), filename.size());
    load_index();
}

void PairedReadsDatastore::load_index(){
    fd=fopen(filename.c_str(),"r");
    fread( &readsize,sizeof(readsize),1,fd);
    fread( &_size,sizeof(_size),1,fd);
    readpos_offset=ftell(fd);
    sglib::OutputLog()<<"PairedReadsDatastore open: "<<filename<<"  max read length: "<<readsize<<" Total reads: " <<size()<<std::endl;
}

void PairedReadsDatastore::load_from_stream(std::string _filename,std::ifstream & input_file){
    uint64_t s;
    filename=_filename;
    fd=fopen(filename.c_str(),"r");
    input_file.read( (char *) &readsize,sizeof(readsize));
    input_file.read( (char *) &_size,sizeof(_size));
    readpos_offset=input_file.tellg();
    fseek(fd,readpos_offset,SEEK_SET);
    input_file.seekg(_size*(readsize+1),std::ios_base::cur);
    sglib::OutputLog()<<"PairedReadsDatastore open: "<<_filename<<"  max read length: "<<readsize<<" Total reads: " <<size()<<std::endl;
}

void PairedReadsDatastore::write(std::ofstream &output_file) {
    //read filename
    uint64_t s=filename.size();
    output_file.write((char *) &s,sizeof(s));
    output_file.write((char *)filename.data(),filename.size());
}

void PairedReadsDatastore::write_selection(std::ofstream &output_file, std::vector<uint64_t> read_ids) {
    for (auto i=0;i<read_ids.size()-1;i+=2){
        if (read_ids[i]+1!=read_ids[i+1]) {
            sglib::OutputLog()<<"ERROR: paired read selection not paired!"<<std::endl;
            return;//exit if not properly paired
        }
    }
    output_file.write((char *) &readsize,sizeof(readsize));
    uint64_t rids_size=read_ids.size();
    output_file.write((char *) &rids_size,sizeof(rids_size));
    char buffer[2*readsize+2];
    for (auto i=0;i<read_ids.size()-1;i+=2) {
        size_t read_offset_in_file = readpos_offset + (readsize + 1) * (read_ids[i] - 1);
        fseek(fd, read_offset_in_file, SEEK_SET);
        fread(buffer, 2 * readsize + 2, 1, fd);
        output_file.write(buffer,2 * readsize + 2);
    }
}

std::string PairedReadsDatastore::get_read_sequence(size_t readID) {
    char buffer[readsize+1];
    size_t read_offset_in_file=readpos_offset+(readsize+1)*(readID-1);
    fseek(fd,read_offset_in_file,SEEK_SET);
    fread(buffer,readsize+1,1,fd);
    return std::string(buffer);
}


const char* BufferedPairedSequenceGetter::get_read_sequence(uint64_t readID) {
    size_t read_offset_in_file=datastore.readpos_offset+(datastore.readsize+1)*(readID-1);
    if (read_offset_in_file<buffer_offset or read_offset_in_file+chunk_size>buffer_offset+bufsize) {
        buffer_offset=read_offset_in_file;
        lseek(fd,read_offset_in_file,SEEK_SET);
        read(fd,buffer,bufsize);
    }
    return buffer+(read_offset_in_file-buffer_offset);
}

namespace std {
    //TODO: this hashing sucks, but it is needed
    template <> struct hash<__int128 unsigned>
    {
        size_t operator()(const __int128 unsigned & x) const
        {
            return hash<uint64_t>()((uint64_t)x);
        }
    };
}

std::unordered_set<__uint128_t> PairedReadsDatastore::get_all_kmers128(int k, int min_tag_cov) {
    class StreamKmerFactory128 : public  KMerFactory128 {
    public:
        explicit StreamKmerFactory128(uint8_t k) : KMerFactory128(k){}
        inline void produce_all_kmers(const char * seq, std::vector<__uint128_t> &mers){
            // TODO: Adjust for when K is larger than what fits in __uint128_t!
            last_unknown=0;
            fkmer=0;
            rkmer=0;
            auto s=seq;
            while (*s!='\0' and *s!='\n') {
                //fkmer: grows from the right (LSB)
                //rkmer: grows from the left (MSB)
                fillKBuf(*s, 0, fkmer, rkmer, last_unknown);
                if (last_unknown >= K) {
                    if (fkmer <= rkmer) {
                        // Is fwd
                        mers.emplace_back(fkmer);
                    } else {
                        // Is bwd
                        mers.emplace_back(rkmer);
                    }
                }
                ++s;
            }
        }
    };
    StreamKmerFactory128 skf(k);

    //reserve space by counting reads first, save only the integer, do not merge just count and insert in the set
    std::vector<__uint128_t> all_kmers;
    BufferedPairedSequenceGetter bprsg(*this,100000,1000);
    for (auto rid=1;rid<=size();++rid) {
        skf.produce_all_kmers(bprsg.get_read_sequence(rid), all_kmers);
    }

    std::sort(all_kmers.begin(),all_kmers.end());
    std::unordered_set<__uint128_t> kset;
    auto ri=all_kmers.begin();
    auto nri=all_kmers.begin();
    while (ri<all_kmers.end()){
        while (nri<all_kmers.end() and *nri==*ri) ++nri;
        if (nri-ri>=min_tag_cov) kset.insert(*ri);
        ri=nri;
    }
    return std::move(kset);
}

std::unordered_set<__uint128_t> PairedReadsDatastore::get_reads_kmers128(int k, int min_tag_cov, std::vector<uint64_t> reads) {
    class StreamKmerFactory128 : public  KMerFactory128 {
    public:
        explicit StreamKmerFactory128(uint8_t k) : KMerFactory128(k){}
        inline void produce_all_kmers(const char * seq, std::vector<__uint128_t> &mers){
            // TODO: Adjust for when K is larger than what fits in __uint128_t!
            last_unknown=0;
            fkmer=0;
            rkmer=0;
            auto s=seq;
            while (*s!='\0' and *s!='\n') {
                //fkmer: grows from the right (LSB)
                //rkmer: grows from the left (MSB)
                fillKBuf(*s, 0, fkmer, rkmer, last_unknown);
                if (last_unknown >= K) {
                    if (fkmer <= rkmer) {
                        // Is fwd
                        mers.emplace_back(fkmer);
                    } else {
                        // Is bwd
                        mers.emplace_back(rkmer);
                    }
                }
                ++s;
            }
        }
    };
    StreamKmerFactory128 skf(k);

    //reserve space by counting reads first, save only the integer, do not merge just count and insert in the set
    std::vector<__uint128_t> all_kmers;
    BufferedPairedSequenceGetter bprsg(*this,100000,1000);
    for (auto rid:reads) {
        skf.produce_all_kmers(bprsg.get_read_sequence(rid), all_kmers);
    }

    std::sort(all_kmers.begin(),all_kmers.end());
    std::unordered_set<__uint128_t> kset;
    auto ri=all_kmers.begin();
    auto nri=all_kmers.begin();
    while (ri<all_kmers.end()){
        while (nri<all_kmers.end() and *nri==*ri) ++nri;
        if (nri-ri>=min_tag_cov) kset.insert(*ri);
        ri=nri;
    }
    return std::move(kset);
}