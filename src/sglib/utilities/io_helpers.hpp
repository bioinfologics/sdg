//
// Created by Bernardo Clavijo (EI) on 28/10/2018.
//

#ifndef BSG_IO_HELPERS_HPP
#define BSG_IO_HELPERS_HPP

namespace sglib {
    template<typename T>
    void write_flat_vector(std::ofstream &ofs, const std::vector<T> &v) {
        uint64_t count=v.size();
        ofs.write(reinterpret_cast<const char *>(&count),sizeof(count));
        ofs.write(reinterpret_cast<const char *>(v.data()),sizeof(T)*count);
    }

    template<typename T>
    void read_flat_vector(std::ifstream &ifs, std::vector<T> &v) {
        uint64_t count;
        ifs.read(reinterpret_cast<char *>(&count),sizeof(count));
        v.resize(count);
        ifs.read(reinterpret_cast<char *>(v.data()),sizeof(T)*count);
    }
}
#endif //BSG_IO_HELPERS_HPP
