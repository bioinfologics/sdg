//
// Created by Bernardo Clavijo (EI) on 2019-06-25.
//

#include "ReadSequenceBuffer.hpp"

ReadSequenceBuffer::ReadSequenceBuffer(const PairedReadsDatastore &_ds, size_t _bufsize , size_t _chunk_size):
        paired_datastore(&_ds),bufsize(_bufsize),chunk_size(_chunk_size){
    fd=open(paired_datastore->filename.c_str(),O_RDONLY);
    buffer=(char *)malloc(bufsize);
    buffer_offset=SIZE_MAX;
}

ReadSequenceBuffer::ReadSequenceBuffer(const LinkedReadsDatastore &_ds, size_t _bufsize , size_t _chunk_size):
        linked_datastore(&_ds),bufsize(_bufsize),chunk_size(_chunk_size){
    fd=open(linked_datastore->filename.c_str(),O_RDONLY);
    buffer=(char *)malloc(bufsize);
    buffer_offset=SIZE_MAX;
}

ReadSequenceBuffer::ReadSequenceBuffer(const LongReadsDatastore &_ds, size_t _bufsize , size_t _chunk_size):
        long_datastore(&_ds),bufsize(_bufsize),chunk_size(_chunk_size){
    fd=open(long_datastore->filename.c_str(),O_RDONLY);
    buffer=(char *)malloc(bufsize);
    buffer_offset=SIZE_MAX;
}

ReadSequenceBuffer::~ReadSequenceBuffer(){
    if (nullptr!=buffer) free(buffer);
    close(fd);
}

const char* ReadSequenceBuffer::get_read_sequence(uint64_t readID) {
    size_t read_offset_in_file;
    if (nullptr!=paired_datastore) {
        read_offset_in_file = paired_datastore->readpos_offset + (paired_datastore->readsize + 1) * (readID - 1);
    } else if (nullptr!=linked_datastore){
        read_offset_in_file=linked_datastore->readpos_offset+(linked_datastore->readsize+1)*(readID-1);
    } else if (nullptr!=long_datastore){
        read_offset_in_file=long_datastore->read_to_fileRecord[readID].offset;
        if (chunk_size < long_datastore->read_to_fileRecord[readID].record_size) {
            throw std::runtime_error(
                    "Reading from " + this->long_datastore->filename +
                    " failed!\nThe size of the buffer chunk is smaller than read " +
                    std::to_string(readID) + " increase the chunk_size so this read fits");
        }
    } else return nullptr;
    if (read_offset_in_file < buffer_offset or read_offset_in_file + chunk_size > buffer_offset + bufsize) {
        buffer_offset = read_offset_in_file;
        lseek(fd, read_offset_in_file, SEEK_SET);
        read(fd, buffer, bufsize);
    }
    return buffer+(read_offset_in_file-buffer_offset);
}