//
// Created by Luis Yanes (EI) on 15/09/2017.
//

#ifndef SEQSORTER_KMERIDX_H
#define SEQSORTER_KMERIDX_H

#include <vector>
#include <limits>
#include <tuple>
#include "sglib/factories/KMerFactory.h"

struct KMerIDXFactoryParams {
    uint8_t k;
};

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
            fillKBuf(currentRecord.seq[p], fkmer, rkmer, last_unknown);
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


class StreamKmerFactory : public  KMerFactory {
public:
    explicit StreamKmerFactory(uint8_t k) : KMerFactory(k){}
    inline void produce_all_kmers(const int64_t seqID, const char * seq, std::vector<KmerIDX> &mers){
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

#endif //SEQSORTER_KMERIDX_H