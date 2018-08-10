#include <iostream>
#include <fstream>
#include <sglib/WorkSpace.hpp>
#include <sglib/processors/GraphMaker.hpp>
#include "sglib/logger/OutputLog.h"
#include "cxxopts.hpp"


int main(int argc, char * argv[]) {
    WorkSpace ws_original;
    std::cout << "Checking " << argv[1] << std::endl;

    ws_original.load_from_disk(argv[1]);

    std::cout << "Num paired datastores " << ws_original.paired_read_datastores.size() << std::endl;
    std::cout << "Num linked datastores " << ws_original.linked_read_datastores.size() << std::endl;
    std::cout << "Num long datastores " << ws_original.long_read_datastores.size() << std::endl;

    std::cout << "Num paired mappers " << ws_original.paired_read_mappers.size() << std::endl;
    std::cout << "Num linked mappers " << ws_original.linked_read_mappers.size() << std::endl;
    std::cout << "Num long mappers " << ws_original.long_read_mappers.size() << std::endl;

    const auto &lrm=ws_original.long_read_mappers[0];
    lrm.update_graph_index();

    uint64_t current_k=lrm.assembly_kmers.begin()->kmer;
    uint32_t count_k(0);
    for (auto it = lrm.assembly_kmers.begin(); it != lrm.assembly_kmers.end(); ++it) {
        if (it->kmer != current_k) {
            std::cout << current_k << "," << count_k << std::endl;
            count_k=0;
            current_k = it->kmer;
        }
        count_k++;
    }

    return 0;
}
