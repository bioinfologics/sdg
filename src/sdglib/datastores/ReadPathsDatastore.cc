//
// Created by Luis Yanes (EI) on 2019-06-26.
//

#include "ReadPathsDatastore.hpp"
#include <fstream>
#include <cstring>
#include <sdglib/logger/OutputLog.hpp>
#include <sdglib/utilities/io_helpers.hpp>

ReadPathsDatastore::ReadPathsDatastore(const std::string &filename) {
    std::ifstream input_file(filename, std::ios::in | std::ios::binary);
    if (!input_file){
        std::cerr << "Couldn't open file: " << filename << ": " << strerror(errno) << std::endl;
        throw std::runtime_error("Couldn't open file: " + filename);
    }
    uint64_t pathcount;
    input_file.read((char *) &pathcount, sizeof(pathcount));
    read_paths.resize(pathcount);
    index.resize(pathcount);
    uint16_t ps;
    uint32_t offset;
    for (uint64_t i = 0; i < pathcount; i++){
        index[i].file_offset = input_file.tellg();
        input_file.read((char *) &offset, sizeof(offset));
        read_paths[i].offset = offset;
        input_file.read((char *) &ps, sizeof(ps));
        sdglib::read_flat_vector(input_file, read_paths[i].path);
    }

    sdglib::OutputLog() << "Loaded ReadPaths file with: " << read_paths.size() << " paths" << std::endl;
}
