//
// Created by Luis Yanes (EI) on 2019-03-21.
//

#ifndef BSG_HASHING_HELPERS_HPP
#define BSG_HASHING_HELPERS_HPP

//  on why XOR is not a good choice for hash-combining:
//  https://stackoverflow.com/questions/5889238/why-is-xor-the-default-way-to-combine-hashes
//
//  this is from boost
//
template<typename T>
inline void hash_combine(std::size_t& seed, const T& val)
{
    std::hash<T> hasher;
    seed ^= hasher(val) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
}

//  taken from https://stackoverflow.com/a/7222201/916549
//
template<typename S, typename T>
struct std::hash<std::pair<S, T>>
{
inline size_t operator()(const std::pair<S, T>& val) const
{
    size_t seed = 0;
    hash_combine(seed, val.first);
    hash_combine(seed, val.second);
    return seed;
}
};

#endif //BSG_HASHING_HELPERS_HPP
