//
// Created by Ben Ward (EI) on 01/06/2020.
//

#include "sdglib/factories/Kmers.h"
#include <string>
#include <iostream>
#include <bitset>

int main() {
    std::cout << "Testing Kmer iteration. Ben J. Ward" << std::endl;
    // A test 75 nucleotide long sequence.

    std::string myseq {"ACCAATCGGGATCGTTTATACCTCTCAGGCGATATTAACTGCGGCGACCCGATTGGAGATACGCTATAGCAAGAT"};

    std::cout << "Testing kmer iteration for string: " << std::endl << myseq << std::endl;

    std::cout << "Creating Kmers<uint64_t,JustKmers>, K = 31" << std::endl;

    Kmers<uint64_t> mykmers(myseq, 31);

    std::cout << "Testing begin iterator creation" << std::endl;
    auto it = mykmers.begin();

    std::cout << "Done" << std::endl << "Testing end iterator creation" << std::endl;
    auto lastit = mykmers.end();

    std::cout << "Done" << std::endl;

    std::cout << "Testing a for loop" << std::endl << std::endl;
    int i {0};
    for(auto it = mykmers.begin(); it != mykmers.end(); it++) {
        std::cout << ++i << ' ' << std::bitset<64>(*it) << std::endl;
    }

    std::cout << std::endl << ".end() scans to find the last iterator position" << std::endl;
    std::cout << "you can also use a condition instead: !it.has_reached_end()" << std::endl << std::endl;

    i = 0;
    for(auto it = mykmers.begin(); !it.has_reached_end(); it++) {
        std::cout << ++i << ' ' << std::bitset<64>(*it) << std::endl;
        std::cout << i << ' ' << std::bitset<64>() << std::endl;
    }

    std::cout << std::bitset<64>(*it) << std::endl;
    std::cout << std::bitset<64>(*lastit) << std::endl;

    std::cout << std::endl << ".end() scans to find the last iterator position" << std::endl;
    std::cout << "you can also use a condition instead: !it.has_reached_end()" << std::endl << std::endl;
}
