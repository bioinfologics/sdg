//
// Created by Ben Ward (EI) on 01/06/2020.
//

#include "sdglib/factories/Kmers.h"
#include <string>
#include <iostream>

int main() {
    // A test 75 nucleotide long sequence.
    std::string myseq {"ACCAATCGGGATCGTTTATACCTCTCAGGCGATATTAACTGCGGCGACCCGATTGGAGATACGCTATAGCAAGAT"};
    Kmers<uint64_t> mykmers(myseq, 31);

}