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
#include <sglib/types/KmerTypes.hpp>

int main(int argc, char **argv){

    auto ref_fds = ::open(argv[1], O_RDONLY);
    if (ref_fds<0) {
        perror( std::string(std::string("Failed opening ") + std::string(argv[1])).c_str());
        exit(1);
    }
    auto asm_fds = ::open(argv[2], O_RDONLY);
    if (asm_fds<0) {
        perror(std::string(std::string("Failed opening ") + std::string(argv[2])).c_str());
        exit(1);
    }

    uint64_t size;
    ::read(ref_fds, &size, sizeof(size));
    ::read(asm_fds, &size, sizeof(size));

    auto bufferSize(10000000u);
    auto count_element_from_ref(0u);
    auto next_element_from_ref ( (KmerIDX*) malloc(bufferSize* sizeof(KmerIDX)));
    auto size_element_from_ref=
            ::read(ref_fds, next_element_from_ref, bufferSize* sizeof(KmerIDX)) / sizeof(KmerIDX);

    auto count_element_from_asm(0u);
    auto next_element_from_asm ( (KmerIDX*) malloc(bufferSize* sizeof(KmerIDX)));
    auto size_element_from_asm=
            ::read(asm_fds, next_element_from_asm, bufferSize* sizeof(KmerIDX)) / sizeof(KmerIDX);

    size_t idxRef(0), idxAsm(0);
    while (count_element_from_ref <= size_element_from_ref and count_element_from_asm <= size_element_from_asm) {
        if (next_element_from_asm[count_element_from_asm].kmer == next_element_from_ref[count_element_from_ref].kmer) {
            std::cout << next_element_from_ref[count_element_from_ref] << "\t" << next_element_from_asm[count_element_from_asm] << std::endl;
            count_element_from_asm++; count_element_from_ref++;

            if (count_element_from_ref==size_element_from_ref){
                size_element_from_ref=
                        ::read(ref_fds, next_element_from_ref, bufferSize* sizeof(KmerIDX)) / sizeof(KmerIDX);
                count_element_from_ref=0;
            }

            if (count_element_from_asm==size_element_from_asm){
                size_element_from_asm=
                        ::read(asm_fds, next_element_from_asm, bufferSize* sizeof(KmerIDX)) / sizeof(KmerIDX);
                count_element_from_asm=0;
            }

        } else if (next_element_from_ref[count_element_from_ref].kmer < next_element_from_asm[count_element_from_asm].kmer){
            count_element_from_ref++;
            if (count_element_from_ref==size_element_from_ref){
                size_element_from_ref=
                        ::read(ref_fds, next_element_from_ref, bufferSize* sizeof(KmerIDX)) / sizeof(KmerIDX);
                count_element_from_ref=0;
            }
        }
        else {
            count_element_from_asm++;
            if (count_element_from_asm==size_element_from_asm){
                size_element_from_asm=
                        ::read(asm_fds, next_element_from_asm, bufferSize* sizeof(KmerIDX)) / sizeof(KmerIDX);
                count_element_from_asm=0;
            }
        }
    }
}