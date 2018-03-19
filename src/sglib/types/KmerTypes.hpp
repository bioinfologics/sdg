//
// Created by Luis Yanes (EI) on 16/03/2018.
//

#ifndef BSG_KMERTYPES_HPP
#define BSG_KMERTYPES_HPP

#include <sglib/SequenceGraph.h>

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


struct graphStrandPos{
    sgNodeID_t node = 0;
    int32_t pos = 0;

    graphStrandPos() {}
    graphStrandPos(sgNodeID_t node, int32_t pos) : node(node), pos(pos) {}

    bool operator==(const graphStrandPos &o) const {
        return (std::tie(node, pos) == std::tie(o.node,o.pos));
    }
};

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
