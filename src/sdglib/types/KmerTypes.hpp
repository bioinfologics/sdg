//
// Created by Luis Yanes (EI) on 16/03/2018.
//

#ifndef BSG_KMERTYPES_HPP
#define BSG_KMERTYPES_HPP

#include <sdglib/types/GenericTypes.hpp>

/**
 * Stores: kmer, node, position and count
 */
struct KmerIDX {

    KmerIDX() : kmer(std::numeric_limits<unsigned long long int>::max()), contigID(0), count(0){}
    explicit KmerIDX(uint64_t kmer) : kmer(kmer), contigID(0), count(0) {}

    KmerIDX(uint64_t _kmer, int32_t _contigID, uint32_t pos = 0, uint8_t _count = 0) : kmer(_kmer), contigID(_contigID),
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

    friend std::ostream& operator<<(std::ostream& os, const KmerIDX& kmer) {
        os << kmer.contigID << "\t" << kmer.pos;
        return os;
    }

    friend class byCtgPos;
    struct byCtgPos {
        bool operator()(const KmerIDX &a, const KmerIDX &b) {
            return std::tie(a.contigID, a.pos) < std::tie(b.contigID,b.pos);
        }
    };


    uint64_t kmer;
    int32_t contigID;
    uint32_t pos;
    uint8_t count;
};

struct KmerIDX_hash {
    size_t operator()(const KmerIDX& k) const {
        return k.kmer;
    }
};

/**
 * This type support k > 31
 */
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
struct KMerIDX128_hash{
    size_t operator()(const KmerIDX128& k) const {
        return (uint64_t) k.kmer;
    }
};

/**
 * Stores a signed node,position pair, where the node encodes the direction
 * The sign on the Node translates to: as seen in node sequence (+) and reverse complement of node sequence(-)
 */
struct graphStrandPos{
    sgNodeID_t node = 0;
    uint32_t pos = 0;

    graphStrandPos() {}
    graphStrandPos(sgNodeID_t node, uint32_t pos) : node(node), pos(pos) {}

    bool operator==(const graphStrandPos &o) const {
        return (std::tie(node, pos) == std::tie(o.node,o.pos));
    }
};

struct KmerCount {

    KmerCount() : kmer(std::numeric_limits<unsigned long long int>::max()), count(0){}
    explicit KmerCount(uint64_t kmer) : kmer(kmer), count(0) {}

    KmerCount(uint64_t _kmer, uint8_t _count) : kmer(_kmer), count(_count) {}

    const bool operator<(const KmerCount& other) const {
        return kmer<other.kmer;
    }

    const bool operator>(const KmerCount &other) const {
        return kmer>other.kmer;
    }

    const bool operator==(const KmerCount &other) const {
        return kmer==other.kmer;
    }

    void merge(const KmerCount &other) {
        uint8_t res = count + other.count;
        res |= -(res < count);

        count = res;
    }

    KmerCount max() {
        return {};
    }

    friend std::ostream& operator<<(std::ostream& os, const KmerCount& kmer) {
        os << kmer.kmer << "\t" << kmer.count;
        return os;
    }

    uint64_t kmer;
    uint8_t count;
    struct KmerCount_hash {
        size_t operator()(const KmerCount& k) const {
            return k.kmer;
        }
    };

};

/**
 * Struct to describe a kmer of the graph, usually employed to describe matches during mapping operations
 * See: LongReadRecruiter::map()
 */
struct kmerPos {
    kmerPos() = default;
    kmerPos(uint64_t kmer, uint32_t contigID, int32_t offset) : kmer(kmer),
                                                                    contigID(contigID),offset(offset) {}

    uint64_t kmer = 0;
    int32_t contigID = 0;
    int32_t offset = 0;

    friend class byKmerContigOffset;

    struct byKmerContigOffset {
        bool operator()(const kmerPos &a, const kmerPos &b) {
            auto a_off(std::abs(a.offset));
            auto b_off(std::abs(b.offset));
//            return std::tie(a.kmer, a.contigID) < std::tie(b.kmer, b.contigID);
            return std::tie(a.kmer, a.contigID, a_off) < std::tie(b.kmer, b.contigID, b_off);
        }
    };

    inline bool operator<(uint64_t const &rhs) const { return kmer < rhs; }
    const bool operator==(const kmerPos &a) const { return std::tie(kmer, contigID, offset) == std::tie(a.kmer, a.contigID, a.offset);}
};

struct kmerPos128 {
    kmerPos128() = default;
    kmerPos128(__uint128_t kmer, uint32_t contigID, int32_t offset) : kmer(kmer),
                                                                contigID(contigID),offset(offset) {}

    __uint128_t kmer = 0;
    int32_t contigID = 0;
    int32_t offset = 0;

    friend class byKmerContigOffset;

    struct byKmerContigOffset {
        bool operator()(const kmerPos128 &a, const kmerPos128 &b) {
            auto a_off(std::abs(a.offset));
            auto b_off(std::abs(b.offset));
//            return std::tie(a.kmer, a.contigID) < std::tie(b.kmer, b.contigID);
            return std::tie(a.kmer, a.contigID, a_off) < std::tie(b.kmer, b.contigID, b_off);
        }
    };

    inline bool operator<(__uint128_t const &rhs) const { return kmer < rhs; }
    const bool operator==(const kmerPos &a) const { return std::tie(kmer, contigID, offset) == std::tie(a.kmer, a.contigID, a.offset);}
};

#endif //BSG_KMERTYPES_HPP
