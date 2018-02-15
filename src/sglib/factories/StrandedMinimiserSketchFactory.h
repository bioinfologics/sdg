//
// Created by Luis Yanes (EI) on 14/02/2018.
//

#ifndef BSG_MINIMISERSKETCHFACTORY_H
#define BSG_MINIMISERSKETCHFACTORY_H

#include <limits>
#include <vector>
#include <algorithm>
#include <set>
#include "KMerFactory.h"

class MinPosIDX {
public:
    uint64_t hash = 0;
    int32_t pos = 0;

    MinPosIDX(uint64_t hash, int32_t pos) : hash(hash), pos(pos) {}
    bool operator<(const MinPosIDX &o) const {
        return hash < o.hash;
    }

    friend std::ostream& operator<<(std::ostream& os, const MinPosIDX& m) {
        os << m.hash;
        return os;
    }
};

/**
 * @brief
 * This class generates a stranded minimiser sketch for a DNA sequence, assumes the sequence contains no N's!!
 */
class StrandedMinimiserSketchFactory : public  KMerFactory {
    uint8_t w;

    inline uint64_t hash(uint64_t key) {
        return key;
        uint64_t mask = (1ULL<<2*K);
        key = (~key + (key << 21)) & mask; // key = (key << 21) - key - 1;
        key = key ^ key >> 24;
        key = ((key + (key << 3)) + (key << 8)) & mask; // key * 265
        key = key ^ key >> 14;
        key = ((key + (key << 2)) + (key << 4)) & mask; // key * 21
        key = key ^ key >> 28;
        key = (key + (key << 31)) & mask;
        return key;
    }

public:
    explicit StrandedMinimiserSketchFactory(uint8_t k, uint8_t w) : KMerFactory(k), w(w) {}
    inline void getMinSketch(const std::string &seq, std::set<MinPosIDX> &sketch){
        // TODO: Adjust for when K is larger than what fits in uint64_t!
        last_unknown=0;
        fkmer=0;
        rkmer=0;
        for (unsigned int nt = 0; nt < K; nt++) {
            fillKBuf(seq[nt], fkmer, rkmer, last_unknown);
        }

        int32_t pos = K;
        for (; pos < seq.length()-K-w+1; ++pos) {
            uint64_t min(std::numeric_limits<uint64_t>::max());
            uint64_t pFkmer(fkmer);
            uint64_t pRkmer(rkmer);
            for (unsigned int j = 0; j < w; j++) {
                fillKBuf(seq[pos + j], fkmer, rkmer, last_unknown);
                if (fkmer == rkmer) continue;
                if (fkmer < rkmer) {
                    // Is fwd
                    min = fkmer < min ? fkmer : min;
                } else {
                    // Is bwd
                    min = rkmer < min ? rkmer : min;
                }
            }
            fkmer = pFkmer;
            rkmer = pRkmer;
            for (unsigned int j = 0; j < w; j++) {
                fillKBuf(seq[pos + j], fkmer, rkmer, last_unknown);
                if (fkmer < rkmer and min == fkmer) {
                    // Is fwd
                    sketch.insert(MinPosIDX(hash(min), pos+j));
                    continue;
                }
                if (rkmer < fkmer and min == rkmer) {
                    // Is bwd
                    sketch.insert(MinPosIDX{hash(min), static_cast<int32_t>(-1*(pos+j))});
                }
            }
        }
    }
};

#endif //BSG_MINIMISERSKETCHFACTORY_H
