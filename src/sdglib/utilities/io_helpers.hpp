//
// Created by Bernardo Clavijo (EI) on 28/10/2018.
//

#pragma once
#include <unordered_set>

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
        if (count > 0) ofs.write(reinterpret_cast<const char *>(v.data()),sizeof(T)*count);
    }

    template<typename T>
    inline void write_flat_vectorvector(std::ofstream &ofs, const std::vector<std::vector<T>> &vv) {
        uint64_t count=vv.size();
        ofs.write(reinterpret_cast<const char *>(&count),sizeof(count));
        if (count > 0) for (auto &v:vv) write_flat_vector(ofs,v);
    }

    template<typename T>
    inline void write_bool_vector(std::ofstream &ofs, const std::vector<T> &v) {
        uint64_t count=v.size();
        ofs.write(reinterpret_cast<const char *>(&count),sizeof(count));
        for(auto i = 0; i < count;) {
            unsigned char aggr = 0;
            for(unsigned char mask = 1; mask > 0 && i < count; ++i, mask <<= 1)
                if(v.at(i)) aggr |= mask;
            ofs.write((const char*)&aggr, sizeof(unsigned char));
        }
    }

    template<typename T1,typename T2>
    inline void write_flat_unorderedmap(std::ofstream &ofs, const std::unordered_map<T1,T2> &u) {
        uint64_t count=u.size();
        ofs.write(reinterpret_cast<const char *>(&count),sizeof(count));
        for (auto e:u) {
            ofs.write(reinterpret_cast<const char *>(&e.first),sizeof(T1));
            ofs.write(reinterpret_cast<const char *>(&e.second),sizeof(T2));
        }
    }

    template<typename T1>
    inline void write_flat_unorderedset(std::ofstream &ofs, const std::unordered_set<T1> &c) {
        uint64_t count=c.size();
        ofs.write(reinterpret_cast<const char *>(&count),sizeof(count));
        for (auto e:c) {
            ofs.write(reinterpret_cast<const char *>(&e),sizeof(T1));
        }
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
        if (count > 0) ifs.read(reinterpret_cast<char *>(v.data()),sizeof(T)*count);
    }

    template<typename T>
    inline void read_flat_vectorvector(std::ifstream &ifs, std::vector<std::vector<T>> &vv) {
        uint64_t count;
        ifs.read(reinterpret_cast<char *>(&count),sizeof(count));
        vv.resize(count);
        if (count > 0) for (auto &v:vv) read_flat_vector<T>(ifs,v);
    }

    template<typename T>
    inline void read_bool_vector(std::ifstream &ifs, std::vector<T> &v) {
        uint64_t count;
        ifs.read(reinterpret_cast<char *>(&count),sizeof(count));
        v.resize(count);
        for(std::vector<bool>::size_type i = 0; i < count;)
        {
            unsigned char aggr;
            ifs.read((char*)&aggr, sizeof(unsigned char));
            for(unsigned char mask = 1; mask > 0 && i < count; ++i, mask <<= 1)
                v.at(i) = aggr & mask;
        }
        //if (count > 0) ifs.read(reinterpret_cast<char *>(v.data()),sizeof(T)*count);
    }

    template<typename T1,typename T2>
    inline void read_flat_unorderedmap(std::ifstream &ifs, std::unordered_map<T1,T2> &u) {
        uint64_t count;
        ifs.read(reinterpret_cast<char *>(&count),sizeof(count));
        u.clear();
        u.reserve(count);
        T1 key;
        T2 value;
        for (auto i=0;i<count;++i) {
            ifs.read(reinterpret_cast<char *>(&key),sizeof(key));
            ifs.read(reinterpret_cast<char *>(&value),sizeof(value));
            u[key]=value;
        }
    }

    template<typename T1>
    inline void read_flat_unorderedset(std::ifstream &ifs, std::unordered_set<T1> &u) {
        uint64_t count;
        ifs.read(reinterpret_cast<char *>(&count),sizeof(count));
        u.clear();
        u.reserve(count);
        T1 value;
        for (auto i=0;i<count;++i) {
            ifs.read(reinterpret_cast<char *>(&value),sizeof(value));
            u.insert(value);
        }
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

