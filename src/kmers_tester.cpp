//
// Created by Ben Ward (EI) on 01/06/2020.
//

#include "sdglib/factories/Kmers.h"
#include "sdglib/factories/KMerFactory.hpp"
#include <string>
#include <iostream>
#include <bitset>
#include <chrono>

int main() {
    std::cout << "Testing Kmer iteration. Ben J. Ward" << std::endl;
    // A test 75 nucleotide long sequence.

    std::string myseq {"ACCAATCGGGATCGTTTATACCTCTCAGGCGATATTAACTGCGGCGACCCGATTGGAGATACGCTATAGCAAGAT"};

    std::cout << "Testing kmer iteration for string: " << std::endl << myseq << std::endl;

    std::cout << "Creating Kmers<uint64_t,JustKmers>, K = 31" << std::endl;

    auto start_stamp = std::chrono::system_clock::now();

    Kmers<uint64_t> mykmers(myseq, 31);

    std::cout << "Testing begin iterator creation" << std::endl;
    auto it = mykmers.begin();

    std::cout << "Done" << std::endl << "Testing end iterator creation" << std::endl;
    auto lastit = mykmers.end();

    std::cout << "Testing a for loop" << std::endl << std::endl;
    for(auto it = mykmers.begin(); it != lastit; it++) {
        std::cout << it->ckmer << std::endl;
    }

    auto end_stamp = std::chrono::system_clock::now();

    std::chrono::duration<double> elapsed_seconds = end_stamp - start_stamp;
    std::time_t end_time = std::chrono::system_clock::to_time_t(end_stamp);

    std::cout << "Finished computation at " << std::ctime(&end_time)
              << "Elapsed time: " << elapsed_seconds.count() << "s\n";

    std::cout << std::endl << ".end() scans to find the last iterator position" << std::endl;
    std::cout << "you can also use a condition instead: !it.has_reached_end()" << std::endl << std::endl;

    auto start_stampb = std::chrono::system_clock::now();


    StringKMerFactory stringfac(31);
    std::vector<uint64_t> mers;
    stringfac.create_kmers(myseq, mers);
    for(auto it = mers.begin(); it != mers.end(); it++) {
        std::cout << *it << std::endl;
    }

    auto end_stampb = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsed_secondsb = end_stampb - start_stampb;
    std::time_t end_timeb = std::chrono::system_clock::to_time_t(end_stampb);

    std::cout << "Finished computation at " << std::ctime(&end_timeb)
              << "Elapsed time: " << elapsed_secondsb.count() << "s\n";

}
