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

/**
 * @brief
 * This class generates a stranded minimiser sketch for a DNA sequence, assumes the sequence contains no N's!!
 */
class StrandedMinimiserSketchFactory : public  KMerFactory {
    uint8_t w;      /// This is the size of the window for minimisers
    uint64_t mask;  /// This masks the bits which are used for the kmers (1<<2*k)-1

    uint64_t hash(uint64_t key) {
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
    explicit StrandedMinimiserSketchFactory(uint8_t k, uint8_t w) : KMerFactory(k), w(w), mask((1ULL<<2*k) - 1) {}

    /**
     * Generate a set of hash(kmer) -> signed position, where the sign represents whether the kmer is forward(+) or reverse(-)
     *
     * @param seq DNA sequence
     * @param sketch Return pre-allocated storage (should be cleared outside this method)
     * @return
     * Returns the number of elements in the "sketch" set.
     */
    inline std::vector<MinPosIDX>::size_type getMinSketch(const std::string &seq, std::vector<MinPosIDX> &sketch){
        // TODO: Adjust for when K is larger than what fits in uint64_t!
        last_unknown=0;
        fkmer=0;
        rkmer=0;
        sketch.reserve(seq.length()/w);
        if (seq.length() < K) return sketch.size();

        for (unsigned int nt = 0; nt < K; nt++) {
            fillKBuf(seq[nt], fkmer, rkmer, last_unknown);
        }

        int32_t pos = K;
        for (; pos < seq.length()-w+1; ++pos) {
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
                    sketch.emplace_back(hash(min), pos+j);
                    continue;
                }
                if (rkmer < fkmer and min == rkmer) {
                    // Is bwd
                    sketch.emplace_back(hash(min), -1*(pos+j));
                }
            }
        }
        return sketch.size();
    }
};

#endif //BSG_MINIMISERSKETCHFACTORY_H
