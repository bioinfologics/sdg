//
// Created by Bernardo Clavijo (EI) on 19/03/2018.
//

#include "PathsDatastore.hpp"

void PathsDatastore::write(std::ofstream &output_file) {
    uint64_t count;
    count=paths.size();
    output_file.write((char *)&count,sizeof(count));
    for (auto &p:paths){
        count=p.nodes.size();
        output_file.write((char *)&count,sizeof(count));
        output_file.write((char *)p.nodes.data(),count*sizeof(p.nodes[0]));
    }
    count=origin.size();
    output_file.write((char *)&count,sizeof(count));
    output_file.write((char *)origin.data(),count*sizeof(origin[0]));
}


void PathsDatastore::read(std::ifstream &input_file) {
    uint64_t count;
    input_file.read((char *)&count,sizeof(count));
    if (input_file.eof()) return;
    SequenceGraphPath p(sg);
    paths.resize(count,p);
    for (auto &p:paths){
        input_file.read((char *)&count,sizeof(count));
        p.nodes.resize(count);
        input_file.read((char *)p.nodes.data(),count*sizeof(p.nodes[0]));
    }
    input_file.read((char *)&count,sizeof(count));
    origin.resize(count);
    input_file.read((char *)origin.data(),count*sizeof(origin[0]));
}