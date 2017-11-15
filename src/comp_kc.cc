//
// Created by Luis Yanes (EI) on 14/11/2017.
//

#include <cstdio>
#include <cstdlib>
#include <set>
#include <unordered_set>
#include <cmath>
#include <fstream>
#include <vector>
#include <limits>
#include <iterator>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>
#include "sglib/KMerIDX.h"


off_t getFilesize(const std::string& filename) {
    struct stat st;
    if(stat(filename.c_str(), &st) != 0) {
        return 0;
    }
    return st.st_size;
}

std::pair<KmerIDX*,uint64_t> readKC(std::string filename) {

    auto in_fds = ::open(filename.c_str(), O_RDONLY);
    uint64_t size;
    ::read(in_fds, &size, sizeof(size));
    // reserve capacity
    KmerIDX* vec = (KmerIDX*) malloc(size*sizeof(KmerIDX));
    // read the data:
    uint64_t numElements = ::read(in_fds, (char *) vec, size*sizeof(KmerIDX)) / sizeof(KmerIDX);
    return std::make_pair(vec,numElements);
}

int main(int argc, char **argv){

    std::pair<KmerIDX*,uint64_t> ref_idx(readKC(std::string(argv[1])));
    std::pair<KmerIDX*,uint64_t> asm_idx(readKC(std::string(argv[2])));

    size_t idxRef(0), idxAsm(0);
    while (idxRef<ref_idx.second and idxAsm < asm_idx.second) {
        if (asm_idx.first[idxAsm].kmer == ref_idx.first[idxRef].kmer) {
            std::cout << ref_idx.first[idxRef] << "\t" << asm_idx.first[idxAsm] << std::endl;
            idxAsm++; idxRef++;
        } else if (ref_idx.first[idxRef].kmer < asm_idx.first[idxAsm].kmer){
            idxRef++;
        }
        else {
            idxAsm++;
        }
    }
}