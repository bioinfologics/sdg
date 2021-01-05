//
// Created by Bernardo Clavijo (EI) on 19/10/2020.
//

#pragma once


#include <sdglib/views/NodeView.hpp>

/**
 * Class to work with populations kmer frequencies
 */
class CountFilter {
public:
    /**
     *
     * @param kcname Name of the KmerCounter to use (needs a computed KmerCounter)
     * @param filter_count_name Name of the count to be used as filter, this restrict if a kmer in the node is
     * considered
     * @param filter_count_max maximum frequency of a kmer in the filter count to be considered
     * @param value_count_names List of names of the counts in the KmerCounter to use in the pattern generation
     * @param value_count_mins Minimum kmer frequency in the counts for a kmer to be considered
     */
    CountFilter(std::string kcname, std::string filter_count_name, int filter_count_max, const std::vector<std::string> value_count_names,const std::vector<int> value_count_mins);

    /** @brief Get the kmer pattern for a node
     *
     * If the frequency of at least kmer in the node is less than value_count_mins for the count being cosidered
     * that node is considered absent from that count and receives a 0 in the pattern, otherwise gets a 1.
     *
     * The results for each count is composed into a string of 1s and 0s following the order of the value_count_name
     * list and returned.
     *
     * @param nv selected node nodeview
     * @return string of the presence/absence of a kmer in the counts in the counts
     */
    std::string get_pattern(NodeView nv);

    std::map<sgNodeID_t,std::string> patterns;
    std::string kcname;
    std::string filter_count_name;
    int filter_count_max;
    std::vector<std::string> value_count_names;
    std::vector<int> value_count_mins;
};



