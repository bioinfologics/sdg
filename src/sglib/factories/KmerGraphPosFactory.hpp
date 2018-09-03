//
// Created by Luis Yanes (EI) on 17/08/2018.
//

#ifndef BSG_KMERGRAPHPOSFACTORY_HPP
#define BSG_KMERGRAPHPOSFACTORY_HPP

#include <vector>
#include "KMerFactory.h"

class KmerGraphPosFactory : protected KMerFactory {
public:
    explicit KmerGraphPosFactory(uint8_t k) : KMerFactory(k) {}

    void setFileRecord(FastaRecord &rec) {
        currentRecord = rec;
        fkmer=0;
        rkmer=0;
        last_unknown=0;
    }

    // TODO: Adjust for when K is larger than what fits in uint64_t!
    const bool next_element(std::vector<std::pair<uint64_t,graphPosition>> &mers) {
        uint64_t p(0);
        graphPosition pos;
        while (p < currentRecord.seq.size()) {
            //fkmer: grows from the right (LSB)
            //rkmer: grows from the left (MSB)
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
};

class KmerGraphPosFactory63 : protected KMerFactory128 {
public:
    explicit KmerGraphPosFactory63(uint8_t k) : KMerFactory128(k) {}

    void setFileRecord(FastaRecord &rec) {
        currentRecord = rec;
        fkmer=0;
        rkmer=0;
        last_unknown=0;
    }

    // TODO: Adjust for when K is larger than what fits in uint64_t!
    const bool next_element(std::vector<std::pair<__uint128_t,graphPosition>> &mers) {
        uint64_t p(0);
        graphPosition pos;
        while (p < currentRecord.seq.size()) {
            //fkmer: grows from the right (LSB)
            //rkmer: grows from the left (MSB)
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
};

#endif //BSG_KMERGRAPHPOSFACTORY_HPP
