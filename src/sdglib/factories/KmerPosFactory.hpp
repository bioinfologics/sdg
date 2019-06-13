//
// Created by Luis Yanes (EI) on 17/08/2018.
//

#ifndef BSG_KMERPOSFACTORY_HPP
#define BSG_KMERPOSFACTORY_HPP

#include <vector>
#include "KMerFactory.hpp"

class kmerPosFactory : protected KMerFactory {
public:
    explicit kmerPosFactory(uint8_t k) : KMerFactory(k) {}

    ~kmerPosFactory() {
    }
    void setFileRecord(FastaRecord &rec) {
        currentRecord = rec;
        fkmer=0;
        rkmer=0;
        last_unknown=0;
    }

    // TODO: Adjust for when K is larger than what fits in uint64_t!
    const bool next_element(std::vector<std::pair<uint64_t,graphStrandPos>> &mers) {
        uint64_t p(0);
        graphStrandPos pos;
        while (p < currentRecord.seq.size()) {
            //fkmer: grows from the right (LSB)
            //rkmer: grows from the left (MSB)
            bases++;
            fillKBuf(currentRecord.seq[p], fkmer, rkmer, last_unknown);
            p++;
            if (last_unknown >= K) {
                if (fkmer <= rkmer) {
                    // Is fwd
                    pos.node=currentRecord.id;
                    pos.pos=p;
                    mers.emplace_back(fkmer, pos);
                } else {
                    // Is bwd
                    pos.node=-currentRecord.id;
                    pos.pos=p;
                    mers.emplace_back(rkmer, pos);
                }
            }
        }
        return false;
    }

private:
    FastaRecord currentRecord;
    uint64_t bases;
};

class kmerPosFactory128 : protected KMerFactory128 {
public:
    explicit kmerPosFactory128(uint8_t k) : KMerFactory128(k) {}

    ~kmerPosFactory128() {
    }
    void setFileRecord(FastaRecord &rec) {
        currentRecord = rec;
        fkmer=0;
        rkmer=0;
        last_unknown=0;
    }

    // TODO: Adjust for when K is larger than what fits in uint64_t!
    const bool next_element(std::vector<std::pair<__uint128_t,graphStrandPos>> &mers) {
        uint64_t p(0);
        graphStrandPos pos;
        while (p < currentRecord.seq.size()) {
            //fkmer: grows from the right (LSB)
            //rkmer: grows from the left (MSB)
            bases++;
            fillKBuf(currentRecord.seq[p], fkmer, rkmer, last_unknown);
            p++;
            if (last_unknown >= K) {
                if (fkmer <= rkmer) {
                    // Is fwd
                    pos.node=currentRecord.id;
                    pos.pos=p;
                    mers.emplace_back(fkmer, pos);
                } else {
                    // Is bwd
                    pos.node=-currentRecord.id;
                    pos.pos=p;
                    mers.emplace_back(rkmer, pos);
                }
            }
        }
        return false;
    }

private:
    FastaRecord currentRecord;
    uint64_t bases;
};

#endif //BSG_KMERPOSFACTORY_HPP
