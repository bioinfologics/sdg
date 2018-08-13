//
// Created by Luis Yanes (EI) on 15/09/2017.
//

#ifndef SEQSORTER_KMERIDX_H
#define SEQSORTER_KMERIDX_H

#include <vector>
#include <limits>
#include <tuple>
#include <cmath>
#include "sglib/factories/KMerFactory.h"
struct KMerIDXFactoryParams {
    uint8_t k;
};
struct kmerPos {
    kmerPos() = default;
    kmerPos(uint64_t kmer, uint32_t contigID, int32_t pos):kmer(kmer),contigID(contigID),pos(pos) {}
    uint64_t  kmer = 0;
    int32_t contigID = 0;
    int32_t pos = 0;

    friend class byKmerContigOffset;
    struct byKmerContigOffset {
        bool operator()(const kmerPos &a, const kmerPos &b) {
            auto a_pos(std::abs(a.pos));
            auto b_pos(std::abs(b.pos));
//            return std::tie(a.kmer, a.contigID) < std::tie(b.kmer, b.contigID);
            return std::tie(a.kmer, a.contigID, a_pos) < std::tie(b.kmer, b.contigID, b_pos);
        }
    };
    inline bool operator<(const uint32_t &km) const {return kmer < km;}
};

struct KmerIDX {

    KmerIDX() : kmer(std::numeric_limits<unsigned long long int>::max()), contigID(0), count(0){}
    explicit KmerIDX(uint64_t kmer) : kmer(kmer), contigID(0), count(0) {}

    KmerIDX(uint64_t _kmer, int32_t _contigID, uint32_t pos, uint8_t _count) : kmer(_kmer), contigID(_contigID),
                                                                               pos(pos), count(_count) { std::cout << contigID << "," << pos << std::endl; }

    const bool operator<(const KmerIDX& other) const {
        return kmer<other.kmer;
    }

    const bool operator>(const KmerIDX &other) const {
        return kmer>other.kmer;
    }

    const bool operator==(const KmerIDX &other) const {
        return kmer==other.kmer;
    }
    void merge(const KmerIDX &other) {
        count += other.count;
    }

    KmerIDX max() {
        return {};
    }

    friend std::ostream& operator<<(std::ostream& os, const KmerIDX& kmer) {
        os << kmer.contigID << "\t" << kmer.pos;
        return os;
    }

    friend std::istream& operator>>(std::istream& is, const KmerIDX& kmer) {
        is.read((char*)&kmer, sizeof(kmer));
        return is;
    }

    friend class byCtgPos;
    struct byCtgPos {
        bool operator()(const KmerIDX &a, const KmerIDX &b) {
            return std::tie(a.contigID, a.pos) < std::tie(b.contigID,b.pos);
        }
    };

    friend class ltKmer;
    struct ltKmer {
        bool operator()(const KmerIDX &a, const uint64_t &b_kmer) const {
            return a.kmer < b_kmer;
        }
    };

    friend class byKmerContigOffset;
    struct byKmerContigOffset {
        bool operator()(const KmerIDX &a, const KmerIDX &b) {
            auto a_pos(std::abs(a.pos));
            auto b_pos(std::abs(b.pos));
//            return std::tie(a.kmer, a.contigID) < std::tie(b.kmer, b.contigID);
            return std::tie(a.kmer, a.contigID, a_pos) < std::tie(b.kmer, b.contigID, b_pos);
        }
    };


    uint64_t kmer;
    int32_t contigID;
    uint32_t pos;
    uint8_t count;
};

struct KmerIDX128 {

    KmerIDX128() : kmer(std::numeric_limits<unsigned long long int>::max()), contigID(0), count(0){}
    explicit KmerIDX128(__uint128_t kmer) : kmer(kmer), contigID(0), count(0) {}

    KmerIDX128(__uint128_t _kmer, int32_t _contigID, uint32_t pos, uint8_t _count) : kmer(_kmer), contigID(_contigID),
                                                                               pos(pos), count(_count) {}

    const bool operator<(const KmerIDX& other) const {
        return kmer<other.kmer;
    }

    const bool operator>(const KmerIDX &other) const {
        return kmer>other.kmer;
    }

    const bool operator==(const KmerIDX &other) const {
        return kmer==other.kmer;
    }
    void merge(const KmerIDX &other) {
        count += other.count;
    }

    KmerIDX max() {
        return {};
    }

    friend std::ostream& operator<<(std::ostream& os, const KmerIDX128& kmer) {
        os << kmer.contigID << "\t" << kmer.pos;
        return os;
    }

    friend std::istream& operator>>(std::istream& is, const KmerIDX128& kmer) {
        is.read((char*)&kmer, sizeof(kmer));
        return is;
    }

    friend class byCtgPos;
    struct byCtgPos {
        bool operator()(const KmerIDX128 &a, const KmerIDX128 &b) {
            return std::tie(a.contigID, a.pos) < std::tie(b.contigID,b.pos);
        }
    };

    __uint128_t kmer;
    int32_t contigID;
    uint32_t pos;
    uint8_t count;
};

namespace std {
    template <>
    struct hash<KmerIDX128> {
        size_t operator()(const KmerIDX128& k) const {
            return (uint64_t) k.kmer;
        }
    };
}

namespace std {
    template <>
    struct hash<KmerIDX> {
        size_t operator()(const KmerIDX& k) const {
            return k.kmer;
        }
    };
}

template<typename FileRecord>
class kmerIDXFactory : protected KMerFactory {
public:
    explicit kmerIDXFactory(KMerIDXFactoryParams params) : KMerFactory(params.k) {}

    ~kmerIDXFactory() {
#pragma omp critical (KMerFactoryDestructor)
        {
            //std::cout << "Bases processed " << bases << "\n";
        }
    }
    void setFileRecord(FileRecord &rec) {
        currentRecord = rec;
        fkmer=0;
        rkmer=0;
        last_unknown=0;
    }

    // TODO: Adjust for when K is larger than what fits in uint64_t!
    const bool next_element(std::vector<KmerIDX> &mers) {
        uint64_t p(0);
        while (p < currentRecord.seq.size()) {
            //fkmer: grows from the right (LSB)
            //rkmer: grows from the left (MSB)
            bases++;
            fillKBuf(currentRecord.seq[p], p, fkmer, rkmer, last_unknown);
            p++;
            if (last_unknown >= K) {
                if (fkmer <= rkmer) {
                    // Is fwd
                    mers.emplace_back(fkmer, currentRecord.id, p, 1);
                } else {
                    // Is bwd
                    mers.emplace_back(rkmer, -1 * currentRecord.id, p, 1);
                }
            }
        }
        return false;
    }

private:
    FileRecord currentRecord;
    uint64_t bases;
};
#endif //SEQSORTER_KMERIDX_H