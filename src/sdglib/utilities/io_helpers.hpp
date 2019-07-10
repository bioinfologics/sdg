//
// Created by Bernardo Clavijo (EI) on 28/10/2018.
//

#pragma once

namespace sdglib {

    inline void write_string(std::fstream &ofs, const std::string &v) {
        uint64_t count=v.size();
        ofs.write(reinterpret_cast<const char *>(&count),sizeof(count));
        ofs.write(reinterpret_cast<const char *>(v.data()),count);
    }

    inline void write_stringvector(std::fstream &ofs, const std::vector<std::string> &vv) {
        uint64_t count=vv.size();
        ofs.write(reinterpret_cast<const char *>(&count),sizeof(count));
        for (auto &v:vv) write_string(ofs,v);
    }

    template<typename T>
    inline void write_flat_vector(std::fstream &ofs, const std::vector<T> &v) {
        uint64_t count=v.size();
        ofs.write(reinterpret_cast<const char *>(&count),sizeof(count));
        ofs.write(reinterpret_cast<const char *>(v.data()),sizeof(T)*count);
    }

    template<typename T>
    inline void write_flat_vectorvector(std::fstream &ofs, const std::vector<std::vector<T>> &vv) {
        uint64_t count=vv.size();
        ofs.write(reinterpret_cast<const char *>(&count),sizeof(count));
        for (auto &v:vv) write_flat_vector(ofs,v);
    }

    inline void write_string(std::ofstream &ofs, const std::string &v) {
        uint64_t count=v.size();
        ofs.write(reinterpret_cast<const char *>(&count),sizeof(count));
        ofs.write(reinterpret_cast<const char *>(v.data()),count);
    }

    inline void write_stringvector(std::ofstream &ofs, const std::vector<std::string> &vv) {
        uint64_t count=vv.size();
        ofs.write(reinterpret_cast<const char *>(&count),sizeof(count));
        for (auto &v:vv) write_string(ofs,v);
    }

    template<typename T>
    inline void write_flat_vector(std::ofstream &ofs, const std::vector<T> &v) {
        uint64_t count=v.size();
        ofs.write(reinterpret_cast<const char *>(&count),sizeof(count));
        ofs.write(reinterpret_cast<const char *>(v.data()),sizeof(T)*count);
    }

    template<typename T>
    inline void write_flat_vectorvector(std::ofstream &ofs, const std::vector<std::vector<T>> &vv) {
        uint64_t count=vv.size();
        ofs.write(reinterpret_cast<const char *>(&count),sizeof(count));
        for (auto &v:vv) write_flat_vector(ofs,v);
    }

    inline void read_string(std::ifstream &ifs, std::string &v) {
        uint64_t count;
        ifs.read(reinterpret_cast<char *>(&count),sizeof(count));
        v.resize(count);
        ifs.read((char *) v.data(),count);
    }

    inline void read_stringvector(std::ifstream &ifs, std::vector<std::string> &vv) {
        uint64_t count;
        ifs.read(reinterpret_cast<char *>(&count),sizeof(count));
        vv.resize(count);
        for (auto &v:vv) read_string(ifs,v);
    }

    template<typename T>
    inline void read_flat_vector(std::ifstream &ifs, std::vector<T> &v) {
        uint64_t count;
        ifs.read(reinterpret_cast<char *>(&count),sizeof(count));
        v.resize(count);
        ifs.read(reinterpret_cast<char *>(v.data()),sizeof(T)*count);
    }

    template<typename T>
    inline void read_flat_vectorvector(std::ifstream &ifs, std::vector<std::vector<T>> &vv) {
        uint64_t count;
        ifs.read(reinterpret_cast<char *>(&count),sizeof(count));
        vv.resize(count);
        for (auto &v:vv) read_flat_vector<T>(ifs,v);
    }

    inline void seek_string(std::fstream &wsfile) {
        uint64_t count2;
        wsfile.read((char *) &count2, sizeof(count2));
        wsfile.seekg(count2, std::ios_base::cur);
    }

    template<typename T>
    inline void seek_flat_vector(std::fstream &ifs, T type) {
        uint64_t count;
        ifs.read(reinterpret_cast<char *>(&count),sizeof(count));
        ifs.seekg(sizeof(type)*count, std::ios_base::cur);
    }

    template<typename T>
    inline void seek_flat_vectorvector(std::fstream &ifs, T type) {
        uint64_t count;
        ifs.read(reinterpret_cast<char *>(&count),sizeof(count));
        for (int i =0; i < count; i++){
            seek_flat_vector<T>(ifs,type);
        }
    }

}

