//
// Created by Luis Yanes (EI) on 18/08/2018.
//

#ifndef BSG_BLOOMFILTER_HPP
#define BSG_BLOOMFILTER_HPP

#include <vector>
#include <cstdint>
#include <cmath>

class BloomFilter {
    uint64_t estimated_size;
    uint32_t numHashes;
    short dataSz = sizeof(data[0])*8;
    std::vector<unsigned char> data;
public:
    BloomFilter(uint64_t estimated_size, uint32_t num_hashes=1) :
    estimated_size(estimated_size),
    numHashes(num_hashes),
    data((estimated_size + 7)/dataSz)
    {}

    void add(const std::vector<uint64_t> &element_hashes) {
        for (const auto& element : element_hashes) {
            auto vector_position(element % estimated_size);
            auto array_position(vector_position / dataSz);
            auto bit_push((7 - vector_position % dataSz));
            auto byte_position_value(1 << bit_push);
#pragma atomic update
            data[array_position] |= byte_position_value;
        }
    }

    void add(const uint64_t &element) {
        auto vector_position(element % estimated_size);
        auto array_position(vector_position / dataSz);
        auto bit_push((7 - vector_position % dataSz));
        auto byte_position_value(1 << bit_push);
#pragma atomic update
        data[array_position] |= byte_position_value;
    }

    bool contains(const std::vector<uint64_t> &element_hashes) const {
        for (const auto& element : element_hashes) {
            auto vector_position(element % estimated_size);
            auto array_position(vector_position / dataSz);
            auto bit_push((7 - vector_position % dataSz));
            auto byte_position_value(1 << bit_push);
            auto result(data[array_position] & byte_position_value);
            if (data[array_position] & byte_position_value != 0)
                return false;
        }
        return true;
    }

    bool contains(const uint64_t &element) const {
        auto vector_position(element % estimated_size);
        auto array_position(vector_position / dataSz);
        auto bit_push((7 - vector_position % dataSz));
        auto byte_position_value(1 << bit_push);
        auto val(data[array_position]);
        auto result(val & byte_position_value);
        return ( result != 0);
    }

    uint64_t number_bits_set() const {
        static const unsigned char setBits[] = {0,1,1,2,1,2,2,3,1,2,2,3,2,3,3,4};
        uint64_t total = 0;
#pragma omp parallel for reduction(+:total)
        for(uint64_t index = 0; index < data.size(); index++) {
            unsigned char current_char = data[index];
            total += setBits[current_char&0x0f]; //Count bottom 4bits
            total += setBits[current_char>>4];  // Count top 4bits
        }
        return total;
    }

    double false_positive_rate() const {
        return pow(double(number_bits_set())/double(estimated_size), numHashes);
    }
};


#endif //BSG_BLOOMFILTER_HPP
