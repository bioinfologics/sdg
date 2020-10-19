//
// Created by Bernardo Clavijo (EI) on 19/10/2020.
//

#pragma once


#include <sdglib/views/NodeView.hpp>

class CountFilter {
public:
    CountFilter(std::string kcname, std::string filter_count_name, int filter_count_max, const std::vector<std::string> value_count_names,const std::vector<int> value_count_mins);
    std::string get_pattern(NodeView nv);
    std::map<sgNodeID_t,std::string> patterns;
    std::string kcname;
    std::string filter_count_name;
    int filter_count_max;
    std::vector<std::string> value_count_names;
    std::vector<int> value_count_mins;
};



