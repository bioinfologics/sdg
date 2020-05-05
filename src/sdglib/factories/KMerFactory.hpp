//
// Created by Luis Yanes (EI) on 20/10/2017.
//

#ifndef SEQ_SORTER_KMERFACTORY_H
#define SEQ_SORTER_KMERFACTORY_H

#include <vector>
#include <cstdint>
#include <iostream>
#include <sdglib/readers/FileReader.hpp>
#include <sdglib/types/KmerTypes.hpp>
#include <array>

#define unlikely(x)     __builtin_expect((x),0)

class KMerFactory {
public:
    const uint8_t K;
    int64_t last_unknown{};
    uint64_t fkmer{};
    uint64_t rkmer{};

    inline void fillKBuf(const unsigned char b, uint64_t &rkmer, uint64_t &fkmer, int64_t &last) {
        if (unlikely(b2f[b] == 4)) {
            last = 0;
            fkmer = ((fkmer << 2) + 0) & KMER_MASK;
            rkmer = (rkmer >> 2) + (((uint64_t) 3) << KMER_FIRSTOFFSET);
        } else {
            fkmer = ((fkmer << 2) + b2f[b]) & KMER_MASK;
            rkmer = (rkmer >> 2) + (((uint64_t) b2r[b]) << KMER_FIRSTOFFSET);
            ++last;
        }
    }

protected:
    explicit KMerFactory(uint8_t k) : K(k), KMER_FIRSTOFFSET((uint64_t) (K - 1) * 2),
                                      KMER_MASK((((uint64_t) 1) << (K * 2)) - 1) {
        for (int i = 0; i < 255; i++) {
            b2f[i] = 4;
            b2r[i] = 4;
        }
        b2f['a'] = b2f['A'] = 0;
        b2f['c'] = b2f['C'] = 1;
        b2f['g'] = b2f['G'] = 2;
        b2f['t'] = b2f['T'] = 3;

        b2r['a'] = b2r['A'] = 3;
        b2r['c'] = b2r['C'] = 2;
        b2r['g'] = b2r['G'] = 1;
        b2r['t'] = b2r['T'] = 0;
    }

private:
    const uint64_t KMER_MASK;
    const uint64_t KMER_FIRSTOFFSET;
    std::array<char,256> b2f={4};
    std::array<char,256> b2r={4};
};

class KMerFactory128 {
public:
    const uint8_t K;
    int64_t last_unknown{};
    __uint128_t fkmer{};
    __uint128_t rkmer{};

    inline void fillKBuf(const unsigned char b, __uint128_t &rkmer, __uint128_t &fkmer, int64_t &last) {
        if (unlikely(b2f[b] == 4)) {
            last = 0;
            fkmer = ((fkmer << 2) + 0) & KMER_MASK;
            rkmer = (rkmer >> 2) + (((__uint128_t) 3) << KMER_FIRSTOFFSET);
        } else {
            fkmer = ((fkmer << 2) + b2f[b]) & KMER_MASK;
            rkmer = (rkmer >> 2) + (((__uint128_t) b2r[b]) << KMER_FIRSTOFFSET);
            ++last;
        }
    }

protected:
    explicit KMerFactory128(uint8_t k) : K(k), KMER_FIRSTOFFSET((__uint128_t) (K - 1) * 2),
                                      KMER_MASK((((__uint128_t) 1) << (K * 2)) - 1) {
        for (int i = 0; i < 255; i++) {
            b2f[i] = 4;
            b2r[i] = 4;
        }
        b2f['a'] = b2f['A'] = 0;
        b2f['c'] = b2f['C'] = 1;
        b2f['g'] = b2f['G'] = 2;
        b2f['t'] = b2f['T'] = 3;

        b2r['a'] = b2r['A'] = 3;
        b2r['c'] = b2r['C'] = 2;
        b2r['g'] = b2r['G'] = 1;
        b2r['t'] = b2r['T'] = 0;
    }

private:
    const __uint128_t KMER_MASK;
    const __uint128_t KMER_FIRSTOFFSET;
    std::array<char,256> b2f={4};
    std::array<char,256> b2r={4};
};


class StringKMerFactory : protected KMerFactory {
public:
    explicit StringKMerFactory(uint8_t k) : KMerFactory(k) {
        fkmer=0;
        rkmer=0;
        last_unknown=0;
    }

    ~StringKMerFactory() {}

    // TODO: Adjust for when K is larger than what fits in uint64_t!
    const bool create_kmers(const std::string &s, std::vector<std::pair<bool, uint64_t>> &mers) {
        fkmer=0;
        rkmer=0;
        last_unknown=0;
        uint64_t p(0);
        mers.reserve(mers.size()+s.size());
        while (p < s.size()) {
            //fkmer: grows from the right (LSB)
            //rkmer: grows from the left (MSB)
            fillKBuf(s[p], fkmer, rkmer, last_unknown);
            p++;
            if (last_unknown >= K) {
                if (fkmer <= rkmer) {
                    // Is fwd
                    mers.emplace_back(true, fkmer);
                } else {
                    // Is bwd
                    mers.emplace_back(false, rkmer);
                }
            }
        }
        return false;
    }

    const bool create_kmers(const std::string &s, std::vector<uint64_t> &mers) {
        fkmer=0;
        rkmer=0;
        last_unknown=0;
        uint64_t p(0);
        mers.reserve(mers.size()+s.size());
        while (p < s.size()) {
            //fkmer: grows from the right (LSB)
            //rkmer: grows from the left (MSB)
            fillKBuf(s[p], fkmer, rkmer, last_unknown);
            p++;
            if (last_unknown >= K) {
                if (fkmer <= rkmer) {
                    // Is fwd
                    mers.emplace_back(fkmer);
                } else {
                    // Is bwd
                    mers.emplace_back(rkmer);
                }
            }
        }
        return false;
    }
};

class StringKMerFactoryNC : protected KMerFactory {
public:
    explicit StringKMerFactoryNC(uint8_t k) : KMerFactory(k) {
        fkmer=0;
        rkmer=0;
        last_unknown=0;
    }

    ~StringKMerFactoryNC() {}

    // TODO: Adjust for when K is larger than what fits in uint64_t!
    const bool create_kmers(const std::string &s, std::vector<std::pair<bool, uint64_t>> &mers) {
        fkmer=0;
        rkmer=0;
        last_unknown=0;
        uint64_t p(0);
        mers.reserve(mers.size()+s.size());
        while (p < s.size()) {
            //fkmer: grows from the right (LSB)
            //rkmer: grows from the left (MSB)
            fillKBuf(s[p], fkmer, rkmer, last_unknown);
            p++;
            if (last_unknown >= K) {
                mers.emplace_back(true, fkmer);
            }
        }
        return false;
    }

    const bool create_kmers(const std::string &s, std::vector<uint64_t> &mers) {
        fkmer=0;
        rkmer=0;
        last_unknown=0;
        uint64_t p(0);
        mers.reserve(mers.size()+s.size());
        while (p < s.size()) {
            //fkmer: grows from the right (LSB)
            //rkmer: grows from the left (MSB)
            fillKBuf(s[p], fkmer, rkmer, last_unknown);
            p++;
            if (last_unknown >= K) {
                mers.emplace_back(fkmer);
            }
        }
        return false;
    }
};


class StreamKmerFactory : protected KMerFactory {
public:
    explicit StreamKmerFactory(uint8_t k) : KMerFactory(k) {}

    inline void produce_all_kmers(const char * seq, std::vector<std::pair<bool, uint64_t>> &mers){
        // TODO: Adjust for when K is larger than what fits in uint64_t!
        last_unknown=0;
        fkmer=0;
        rkmer=0;
        auto s=seq;
        while (*s!='\0' and *s!='\n') {
            //fkmer: grows from the right (LSB)
            //rkmer: grows from the left (MSB)
            fillKBuf(*s, fkmer, rkmer, last_unknown);
            if (last_unknown >= K) {
                if (fkmer <= rkmer) {
                    // Is fwd
                    mers.emplace_back(true,fkmer);
                } else {
                    // Is bwd
                    mers.emplace_back(false,rkmer);
                }
            }
            ++s;
        }
    }

    inline void produce_all_kmers(const char * seq, std::vector<uint64_t> &mers){
        // TODO: Adjust for when K is larger than what fits in uint64_t!
        last_unknown=0;
        fkmer=0;
        rkmer=0;
        auto s=seq;
        while (*s!='\0' and *s!='\n') {
            //fkmer: grows from the right (LSB)
            //rkmer: grows from the left (MSB)
            fillKBuf(*s, fkmer, rkmer, last_unknown);
            if (last_unknown >= K) {
                if (fkmer <= rkmer) {
                    // Is fwd
                    mers.emplace_back(fkmer);
                } else {
                    // Is bwd
                    mers.emplace_back(rkmer);
                }
            }
            ++s;
        }
    }

};

class StreamKmerFactory128 : public  KMerFactory128 {
public:
    explicit StreamKmerFactory128(uint8_t k) : KMerFactory128(k){}
    inline void produce_all_kmers(const char * seq, std::vector<std::pair<bool, __uint128_t>> &mers){
        // TODO: Adjust for when K is larger than what fits in uint64_t!
        last_unknown=0;
        fkmer=0;
        rkmer=0;
        auto s=seq;
        while (*s!='\0' and *s!='\n') {
            //fkmer: grows from the right (LSB)
            //rkmer: grows from the left (MSB)
            fillKBuf(*s, fkmer, rkmer, last_unknown);
            if (last_unknown >= K) {
                if (fkmer <= rkmer) {
                    // Is fwd
                    mers.emplace_back(true,fkmer);
                } else {
                    // Is bwd
                    mers.emplace_back(false,rkmer);
                }
            }
            ++s;
        }
    }
    inline void produce_all_kmers(const char * seq, std::vector<__uint128_t> &mers){
        // TODO: Adjust for when K is larger than what fits in __uint128_t!
        last_unknown=0;
        fkmer=0;
        rkmer=0;
        auto s=seq;
        while (*s!='\0' and *s!='\n') {
            //fkmer: grows from the right (LSB)
            //rkmer: grows from the left (MSB)
            fillKBuf(*s, fkmer, rkmer, last_unknown);
            if (last_unknown >= K) {
                if (fkmer <= rkmer) {
                    // Is fwd
                    mers.emplace_back(fkmer);
                } else {
                    // Is bwd
                    mers.emplace_back(rkmer);
                }
            }
            ++s;
        }
    }
};

class CStringKMerFactory : protected KMerFactory {
public:
    explicit CStringKMerFactory(uint8_t k) : KMerFactory(k) {};



    ~CStringKMerFactory() {
#pragma omp critical (KMerFactoryDestructor)
        {
            //std::cout << "Bases processed " << bases << "\n";
        }
    }

    // TODO: Adjust for when K is larger than what fits in uint64_t!
    const void create_kmercounts(std::vector<KmerCount> &mers, const char * s) {
        fkmer=0;
        rkmer=0;
        last_unknown=0;
        uint64_t p(0);
        while (s[p]!='\0') {
            //fkmer: grows from the right (LSB)
            //rkmer: grows from the left (MSB)
            fillKBuf(s[p], fkmer, rkmer, last_unknown);
            p++;
            if (last_unknown >= K) {
                if (fkmer <= rkmer) {
                    // Is fwd
                    mers.emplace_back(fkmer,1);
                } else {
                    // Is bwd
                    mers.emplace_back(rkmer,1);
                }
            }
        }
    }
    const void create_kmers(std::vector<uint64_t> &mers, const char * s) {
        fkmer=0;
        rkmer=0;
        last_unknown=0;
        uint64_t p(0);
        while (s[p]!='\0') {
            //fkmer: grows from the right (LSB)
            //rkmer: grows from the left (MSB)
            fillKBuf(s[p], fkmer, rkmer, last_unknown);
            p++;
            if (last_unknown >= K) {
                if (fkmer <= rkmer) {
                    // Is fwd
                    mers.emplace_back(fkmer);
                } else {
                    // Is bwd
                    mers.emplace_back(rkmer);
                }
            }
        }
    }
    const void create_kmers_direction(std::vector<std::pair<uint64_t,bool>> &mers, const char * s) {
        fkmer=0;
        rkmer=0;
        last_unknown=0;
        uint64_t p(0);
        while (s[p]!='\0') {
            //fkmer: grows from the right (LSB)
            //rkmer: grows from the left (MSB)
            fillKBuf(s[p], fkmer, rkmer, last_unknown);
            p++;
            if (last_unknown >= K) {
                if (fkmer <= rkmer) {
                    // Is fwd
                    mers.emplace_back(fkmer,false);
                } else {
                    // Is bwd
                    mers.emplace_back(rkmer,true);
                }
            }
        }
    }
};

class CStringKMerFactory128 : protected KMerFactory128 {
public:
    explicit CStringKMerFactory128(uint8_t k) : KMerFactory128(k) {};



    ~CStringKMerFactory128() {
#pragma omp critical (KMerFactoryDestructor)
        {
            //std::cout << "Bases processed " << bases << "\n";
        }
    }

    const void create_kmers_direction(std::vector<std::pair<__uint128_t,bool>> &mers, const char * s) {
        fkmer=0;
        rkmer=0;
        last_unknown=0;
        uint64_t p(0);
        while (s[p]!='\0') {
            //fkmer: grows from the right (LSB)
            //rkmer: grows from the left (MSB)
            fillKBuf(s[p], fkmer, rkmer, last_unknown);
            p++;
            if (last_unknown >= K) {
                if (fkmer <= rkmer) {
                    // Is fwd
                    mers.emplace_back(fkmer,false);
                } else {
                    // Is bwd
                    mers.emplace_back(rkmer,true);
                }
            }
        }
    }
};

namespace sdglib {
    static std::vector<std::pair<bool, uint64_t>> str_to_kmers(std::string s, uint8_t k) {
        std::vector<std::pair<bool, uint64_t>> r;
        StreamKmerFactory skf(k);
        skf.produce_all_kmers(s.c_str(), r);
        return r;
    }
}

#endif //SEQ_SORTER_KMERFACTORY_H
