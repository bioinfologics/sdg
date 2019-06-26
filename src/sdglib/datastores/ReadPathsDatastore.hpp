//
// Created by Luis Yanes (EI) on 2019-06-26.
//

#pragma once

#include <sdglib/types/GenericTypes.hpp>
#include <vector>

struct ReadPath{
    uint32_t offset;
    std::vector<sgNodeID_t> path;
};

struct ReadPathOffset {
    uint64_t file_offset;
};

class ReadPathsDatastore {

public:
    explicit ReadPathsDatastore(const std::string &filename);

    std::vector<ReadPathOffset> index; // Stores the file location for the path object
    std::vector<ReadPath> read_paths;
};