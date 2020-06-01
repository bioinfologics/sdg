//
// Created by Ben Ward (EI) on 01/06/2020.
//

#include "sdglib/factories/Kmers.h"
#include <string>
#include <iostream>

int main() {
    // A test 75 nucleotide long sequence.
    std::string myseq {"ACCAATCGGGATCGTTTATACCTCTCAGGCGATATTAACTGCGGCGACCCGATTGGAGATACGCTATAGCAAGAT"};
    Kmers<uint64_t,JustKmers> mykmers(myseq, 31);

    for(auto it = mykmers.begin(); it < mykmers.end(); ++it) {
        std::cout << *it << std::endl;
    }
}