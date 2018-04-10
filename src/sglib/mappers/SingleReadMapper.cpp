//
// Created by Luis Yanes (EI) on 09/04/2018.
//

#include "SingleReadMapper.hpp"
#include <sglib/utilities/omp_safe.hpp>


void SingleReadMapper::map_reads(const std::unordered_set<uint64_t> &readIDs) {

    std::vector< std::vector<LongReadMapping> > read_mappings(datastore.size());
    if (not readIDs.empty())
        sglib::OutputLog()<<readIDs.size()<<" selected reads / "<<datastore.size()-1<<" total"<<std::endl;

#pragma omp parallel
    {
        std::ofstream matchOutput(std::string("thread_")+std::to_string(omp_get_thread_num())+std::string(".match"));
        std::ofstream blockOutput(std::string("thread_")+std::to_string(omp_get_thread_num())+std::string(".block"));
#pragma omp for
        for (uint32_t readID=1;readID<datastore.size();++readID) {
            if ( (!readIDs.empty() and readIDs.count(readID)) or (readIDs.empty() and read_to_mappings[readID].empty()) ) {
                const auto currentReadMappings(map_read
                                                       (readID,
                                                        datastore.get_read_sequence(readID),
                                                        matchOutput,
                                                        blockOutput)
                );
                read_mappings[readID] = currentReadMappings;
            }
        }
    }
    // Initialise reads_in_node, read_to_nodes from results
    std::vector<LongReadMapping>::size_type totalMappings(0);
    for (const auto &rm : read_mappings) {
        totalMappings += rm.size();
    }
    mappings.reserve(totalMappings);
    for (const auto &rm : read_mappings) {
        mappings.insert(mappings.end(), rm.cbegin(), rm.cend());
    }

    for (std::vector<LongReadMapping>::size_type i=1; i < mappings.size(); ++i) {
        mappings_in_node[std::abs(mappings[i].node)].push_back(i);
        read_to_mappings[mappings[i].read_id].push_back(i);
    }

}
