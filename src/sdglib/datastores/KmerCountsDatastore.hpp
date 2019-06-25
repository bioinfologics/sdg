//
// Created by Bernardo Clavijo (EI) on 2019-06-25.
//

#pragma once


#include <sdglib/workspace/WorkSpace.hpp>

class KmerCountsDatastore {
    KmerCountsDatastore(const WorkSpace &_ws):ws(_ws){};
    void add_count(const std::string & count_name, const std::vector<std::string> & filenames);
    template<class T>
    void add_count(const std::string & count_name, const T & datastore);

    std::vector<std::string> count_names;
    std::vector<std::vector<uint16_t>> counts;
    const WorkSpace &ws;
};
