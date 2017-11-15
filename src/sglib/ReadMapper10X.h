//
// Created by Katie Barr (EI) on 15/11/2017.
//

#ifndef SG_READMAPPER10X_H
#define SG_READMAPPER10X_H


#include <map>
#include<atomic>
#include "SequenceGraph.hpp"
#include "KMerIDX.h"
#include "FileReader.h"
// very similar to PairedReadMapper, but grouped by barcode- call from paired read mapper function

class ReadMapping10X {
public:
    sgNodeID_t node;
    std::string barcode;
    uint32_t read_id;
    int32_t first_pos;
    int32_t last_pos;
    int32_t unique_matches;
    bool rev;
    bool operator==(const ReadMapping10X &other){
        return this==&other;
    };
    bool operator<(const ReadMapping10X &other) const {
        return node<other.node;
    };
    void merge(const ReadMapping10X &other){};
};

class Barcode {
public:
    std::string barcode;
    std::vector<ReadMapping10X> mappings;
};


class ReadMapper10X {
public:
    std::vector<Barcode> barcodes;
    uint64_t process_reads_from_file(uint8_t, uint16_t , std::vector<KmerIDX> &, std::string , uint64_t  );

private:
    std::map<std::string, int> barcode_index_map;

};


#endif //SG_READMAPPER10X_H
