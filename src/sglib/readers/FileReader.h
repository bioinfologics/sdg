//
// Created by Luis Yanes (EI) on 15/09/2017.
//

#ifndef SEQSORTER_FILEREADER_H
#define SEQSORTER_FILEREADER_H

#include <cstring>
#include <fstream>
#include <iostream>
#include <fcntl.h>
#include "Common.h"
#include "kseq.hpp"

struct FastxReaderParams {
    uint32_t min_length;
};

struct FastaRecord{
    int32_t id;
    std::string name,seq;
};

template<typename FileRecord>
class FastaReader {
public:
    /**
     * @brief
     * Initialises the FastaReader, opens the file based on the format and instantiates a reader (plain, gzip or bzip2)
     * @param params
     * Parameters for filtering the records (i.e. min_size, max_size)
     * @param filepath
     * Relative or absolute path to the file that is going to be read.
     */
    explicit FastaReader(FastxReaderParams params, const std::string &filepath) : params(params), numRecords(0) {
        std::cout << "Opening: " << filepath << "\n";
        gz_file = gzopen(filepath.c_str(), "r");
        if (gz_file == Z_NULL) {
            std::cout << "Error opening FASTA " << filepath << ": " << std::strerror(errno) << std::endl;
            exit(1);
        }
        ks = new kstream<gzFile, FunctorZlib>(gz_file, gzr);
    }

    /**
     * @brief
     * Calls the file reader and places the fields from the file onto the FileRecord, the ID is set to the
     * number of records seen so far.
     * @param rec
     * Input/Output parameter where the file fields will be stored.
     * @return
     * Whether the function will generate another object or not
     */
    bool next_record(FileRecord& rec) {
        int l;
        do {
            l=(ks->readFasta(seq));
            std::swap(rec.seq, seq.seq);
            std::swap(rec.name, seq.name);
            rec.id = numRecords;
            numRecords++;
            stats.totalLength+=rec.seq.size();
        } while(rec.seq.size() < params.min_length && l >= 0);
        stats.filteredRecords++;
        stats.filteredLength+=rec.seq.size();
        return (l >= 0);
    }
    ReaderStats getSummaryStatistics() {
        stats.totalRecords = numRecords;
        return stats;
    }
private:
    kstream<gzFile, FunctorZlib> *ks;
    kstream<BZFILE, FunctorBZlib2> *bzKS;
    kseq seq;
    uint32_t numRecords=0;
    gzFile gz_file;
    BZFILE * bz_File{};
    int fq_File{};
    FunctorZlib gzr;
    FunctorRead rr;
    FastxReaderParams params;
    ReaderStats stats;
};

struct FastqRecord{
    int64_t id;
    std::string name, comment, seq, qual;
};

template<typename FileRecord>
class FastqReader {
public:

    /**
     * @brief
     * Initialises the FastaReader, opens the file based on the format and instantiates a reader (plain, gzip or bzip2)
     * @param params
     * Parameters for filtering the records (i.e. min_size, max_size)
     * @param filepath
     * Relative or absolute path to the file that is going to be read.
     */
    explicit FastqReader(FastxReaderParams params, const std::string &filepath) : params(params), numRecords(0),eof_flag(false) {
        std::cout << "Opening: " << filepath << "\n";
        gz_file = gzopen(filepath.c_str(), "r");
        if (gz_file == Z_NULL) {
            std::cout << "Error opening FASTQ " << filepath << ": " << std::strerror(errno) << std::endl;
            exit(1);
        }
        ks = new kstream<gzFile, FunctorZlib>(gz_file, gzr);
    }

    /**
    * @brief
    * Calls the file reader and places the fields from the file onto the FileRecord, the ID is set to the
    * number of records seen so far.
    * @param rec
    * Input/Output parameter where the file fields will be stored.
    * @return
    * Whether the function will generate another object or not
    */
    bool next_record(FileRecord& rec) {
        int l;
        if ( eof_flag) return false;
        {
            do {
                l = (ks->readFastq(seq));
                std::swap(rec.seq, seq.seq);
                std::swap(rec.qual, seq.qual);
                std::swap(rec.name, seq.name);
                std::swap(rec.comment, seq.comment);
                rec.id = numRecords;
                numRecords++;

                stats.totalLength += rec.seq.size();
            } while (rec.seq.size() < params.min_length && l >= 0);
        }
        if (l<0) eof_flag=true;
        else {
            stats.filteredRecords++;
            stats.filteredLength += rec.seq.size();
        }

        return (l >= 0);
    }

    ReaderStats getSummaryStatistics() {
        stats.totalRecords=numRecords;
        return stats;
    }

private:
    kstream<gzFile, FunctorZlib> *ks;
    kseq seq;
    uint64_t numRecords=0;
    gzFile gz_file;
    BZFILE * bz_File{};
    int fq_File{};
    FunctorZlib gzr;
    FastxReaderParams params;
    ReaderStats stats;
    bool eof_flag;
};

#endif //SEQSORTER_FILEREADER_H
