//
// Created by Luis Yanes (EI) on 16/03/2018.
//

#ifndef BSG_KMERTYPES_HPP
#define BSG_KMERTYPES_HPP

#include <sglib/types/GenericTypes.hpp>

/**
 * An SMR compatible kmer type provides: count(field), merge(func) and max(constructor)
 * Stores: kmer, node, position and count
 */
struct KmerIDX {

    KmerIDX() : kmer(std::numeric_limits<unsigned long long int>::max()), contigID(0), count(0){}
    explicit KmerIDX(uint64_t kmer) : kmer(kmer), contigID(0), count(0) {}

    KmerIDX(uint64_t _kmer, int32_t _contigID, uint32_t pos, uint8_t _count) : kmer(_kmer), contigID(_contigID),
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

    uint64_t kmer;
    int32_t contigID;
    uint32_t pos;
    uint8_t count;
};

namespace std {
    template <>
    struct hash<KmerIDX> {
        size_t operator()(const KmerIDX& k) const {
            return k.kmer;
        }
    };
}


/**
 * Stores a ±node,position pair, where the node encodes the direction
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

/**
 * Stores a hash(kmer),±position pair where the position encodes for fwd or rc of the kmer
 * The sign on the position translates to: as seen on the input sequence (+) and reverse complement of the input seq (-)
 */
class MinPosIDX {
public:
    uint64_t hash = 0;
    int32_t pos = 0;

    MinPosIDX() : hash(0), pos(0) {}

    MinPosIDX(uint64_t hash, int32_t pos) : hash(hash), pos(pos) {}
    bool operator<(const MinPosIDX &o) const {
        return hash < o.hash;
    }

    friend std::ostream& operator<<(std::ostream& os, const MinPosIDX& m) {
        os << m.hash;
        return os;
    }
};


struct graphPosition{
    sgNodeID_t node;
    uint32_t pos;
};

#endif //BSG_KMERTYPES_HPP
