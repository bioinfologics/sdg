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
#include "sglib/KMerIDX.h"

struct KmerIDXMatch{
    KmerIDX a,b;
};

std::vector<KmerIDX> readKC(std::string filename) {
    std::ifstream ref_kc(filename, std::ios_base::binary);
    ref_kc.unsetf(std::ios_base::skipws);

    std::streampos fileSize;
    ref_kc.seekg(0, std::ios_base::end);
    fileSize = ref_kc.tellg();
    ref_kc.seekg(sizeof(uint64_t), std::ios_base::beg);
    // reserve capacity
    std::vector<KmerIDX> vec;
    vec.reserve(fileSize);

    // read the data:
    vec.insert(vec.begin(), std::istream_iterator<KmerIDX>(ref_kc), std::istream_iterator<KmerIDX>());
    return vec;
}

int main(int argc, char **argv){

    std::vector<KmerIDX> ref_idx(readKC(std::string(argv[1])));
    std::vector<KmerIDX> asm_idx(readKC(std::string(argv[2])));

    size_t idxRef(0), idxAsm(0);
    while (idxRef<ref_idx.size() and idxAsm < asm_idx.size()) {
        if (asm_idx[idxAsm].kmer == ref_idx[idxRef].kmer) {
            std::cout << ref_idx[idxRef] << "\t" << asm_idx[idxAsm] << std::endl;
            idxAsm++; idxRef++;
        } else if (ref_idx[idxRef].kmer < asm_idx[idxAsm].kmer){
            idxRef++;
        }
        else {
            idxAsm++;
        }
    }
}