//
// Created by Bernardo Clavijo (EI) on 2019-06-25.
//

#pragma once

#include <cstdint>
#include <cstdio>

class PairedReadsDatastore;
class LinkedReadsDatastore;
class LongReadsDatastore;

/**
 * This class accesses the sequence of the reads of a datastore through an internal buffer.
 * It can be constructed with any kind of read datastore and will return a pointer a its internal buffer containing the
 * read's sequence.
 */


class ReadSequenceBuffer {
public:
    explicit ReadSequenceBuffer(const PairedReadsDatastore &_ds, size_t _bufsize = (1024*1024*30ul), size_t _chunk_size = (1024*1024*4ul));
    explicit ReadSequenceBuffer(const LongReadsDatastore &_ds, size_t _bufsize = (1024*1024*30ul), size_t _chunk_size = (1024*1024*4ul));
    explicit ReadSequenceBuffer(const LinkedReadsDatastore &_ds, size_t _bufsize = (1024*1024*30ul), size_t _chunk_size = (1024*1024*4ul));
    const char * get_read_sequence(uint64_t readID);
    ~ReadSequenceBuffer();
    ReadSequenceBuffer& operator=(const ReadSequenceBuffer&) = delete;
private:
    const PairedReadsDatastore * paired_datastore= nullptr;
    const LinkedReadsDatastore * linked_datastore= nullptr;
    const LongReadsDatastore * long_datastore= nullptr;
    char * buffer= nullptr;
    size_t bufsize,chunk_size;
    size_t buffer_offset;
    int fd;
};