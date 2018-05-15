//
// Created by Bernardo Clavijo (EI) on 10/02/2018.
//

#include <sglib/logger/OutputLog.h>
#include <fstream>
#include <strings.h>
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
    uint64_t rts=pairs*2;
    output.write((const char *) &rts, sizeof(rts));
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
        if (currrent_read.seq1.size()>_rs) ++truncated;
        if (currrent_read.seq2.size()>_rs) ++truncated;
        ++pairs;
        readdatav.push_back(currrent_read);
        if (readdatav.size()==chunksize){
            //dump
            char buffer[2*readsize+2];
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
        for (auto &r:readdatav){
            bzero(buffer,2*readsize+2);
            memcpy(buffer,r.seq1.data(),(r.seq1.size()>readsize ? readsize : r.seq1.size()));
            memcpy(buffer+readsize+1,r.seq2.data(),(r.seq2.size()>readsize ? readsize : r.seq2.size()));
            output.write(buffer,2*readsize+2);
        }
        sglib::OutputLog()<<readdatav.size()<<" pairs dumped..."<<std::endl;
        readdatav.clear();
    }
    rts=pairs*2;
    output.seekp(size_pos);
    output.write((const char *) &rts, sizeof(rts));
    output.close();
    //DONE!
    sglib::OutputLog(sglib::LogLevels::INFO)<<discarded<<" pairs discarded due to short reads"<<std::endl;
    sglib::OutputLog(sglib::LogLevels::INFO)<<truncated<<" reads where truncated to "<<_rs<<"bp"<<std::endl;
    sglib::OutputLog(sglib::LogLevels::INFO)<<"Datastore with "<<rts<<" reads ("<<pairs<<" pairs)"<<std::endl;
    filename=output_filename;
    fd=fopen(filename.c_str(),"r");
    _size=rts;
}

void PairedReadsDatastore::read(std::ifstream &input_file) {
    //read filename
    uint64_t s;
    input_file.read((char *) &s, sizeof(s));
    filename.resize(s);
    input_file.read((char *) filename.data(), filename.size());
    load_index(filename);
}

void PairedReadsDatastore::load_index(std::string _filename){
    uint64_t s;
    filename=_filename;
    fd=fopen(filename.c_str(),"r");
    fread( &readsize,sizeof(readsize),1,fd);
    fread(&s,sizeof(s),1,fd);
    _size=s;
    readpos_offset=ftell(fd);
}

void PairedReadsDatastore::write(std::ofstream &output_file) {
    //read filename
    uint64_t s=filename.size();
    output_file.write((char *) &s,sizeof(s));
    output_file.write((char *)filename.data(),filename.size());
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
