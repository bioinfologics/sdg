//
// Created by Bernardo Clavijo (EI) on 19/03/2018.
//

#include <fstream>
#include <sglib/graph/SequenceGraphPath.hpp>
#include "PathsDatastore.hpp"


const bsgVersion_t PathsDatastore::min_compat = 0x0001;
void PathsDatastore::write(std::ofstream &output_file) {
    uint64_t count;
    count=paths.size();
    output_file.write((char *) &BSG_MAGIC, sizeof(BSG_MAGIC));
    output_file.write((char *) &BSG_VN, sizeof(BSG_VN));
    BSG_FILETYPE type(PathDS_FT);
    output_file.write((char *) &type, sizeof(type));

    output_file.write((char *)&count,sizeof(count));
    for (auto &p:paths){
        count=p.getNodes().size();
        output_file.write((char *)&count,sizeof(count));
        output_file.write((char *)p.getNodes().data(),count*sizeof(p.getNodes()[0]));
    }
    count=origin.size();
    output_file.write((char *)&count,sizeof(count));
    output_file.write((char *)origin.data(),count*sizeof(origin[0]));
}


void PathsDatastore::read(std::ifstream &input_file) {
    bsgMagic_t magic;
    bsgVersion_t version;
    BSG_FILETYPE type;
    input_file.read((char *) &magic, sizeof(magic));
    input_file.read((char *) &version, sizeof(version));
    input_file.read((char *) &type, sizeof(type));

    if (magic != BSG_MAGIC) {
        throw std::runtime_error("PathsDatastore file seems to be corrupted");
    }

    if (version < min_compat) {
        throw std::runtime_error("Incompatible version");
    }

    if (type != PathDS_FT) {
        throw std::runtime_error("Incompatible file type");
    }

    uint64_t count;
    input_file.read((char *)&count,sizeof(count));
    if (input_file.eof()) return;
    SequenceGraphPath p(sg);
    paths.resize(count,p);
    for (auto &p:paths){
        input_file.read((char *)&count,sizeof(count));
        p.getNodes().resize(count);
        input_file.read((char *)p.getNodes().data(),count*sizeof(p.getNodes()[0]));
    }
    input_file.read((char *)&count,sizeof(count));
    origin.resize(count);
    input_file.read((char *)origin.data(),count*sizeof(origin[0]));
}

PathsDatastore &PathsDatastore::operator=(const PathsDatastore &other) {
    if (&sg != &other.sg) { throw ("Can only copy paths from the same SequenceGraph"); }
    if (&other == this) {
        return *this;
    }
    origin = other.origin;
    paths = other.paths;
    return *this;
}
