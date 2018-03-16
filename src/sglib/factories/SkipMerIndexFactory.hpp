//
// Created by Luis Yanes (EI) on 16/03/2018.
//

#ifndef BSG_SKIPMERINDEXFACTORY_HPP
#define BSG_SKIPMERINDEXFACTORY_HPP

#include <vector>
#include <sglib/types/KmerTypes.hpp>
#include "SipMerFactory.h"

template<typename FileRecord>
class SkipMerIDXFactory : protected SkipMerFactory {
public:
    explicit SkipMerIDXFactory(uint8_t k, uint8_t m, uint8_t n) : SkipMerFactory(k, m, n) {}

    void setFileRecord(FileRecord &rec) {
        currentRecord = rec;
        for (auto ni=0; ni < N; ++ni) {
            fkmer[ni] =0;
            rkmer[ni] = 0;
            last_unknown[ni] = 0;
        }
        fi=0;
    }

    const bool next_element(FileRecord &currentRecord, std::vector<MinPosIDX> &mers){
        for (auto ni=0; ni < N; ++ni) {
            fkmer[ni] =0;
            rkmer[ni] = 0;
            last_unknown[ni] = 0;
        }
        fi=0;
        uint64_t p(0);
        while (p < currentRecord.seq.size()) {
            //fkmer: grows from the right (LSB)
            //rkmer: grows from the left (MSB)
            bases++;
            for (auto ni=0;ni<N;++ni) {
                cycle_pos[ni]++;
                if (cycle_pos[ni]==N)cycle_pos[ni]=0;
                if (cycle_pos[ni]<M) {
                    fillKBuf(currentRecord.seq[p], fkmer[ni], rkmer[ni], last_unknown[ni]);
                }
                p++;
                //if we are at p, the skip-mer started at p-S is now done
                if (p>=S-1) {
                    if (p == S - 1) fi = 0;
                    else {
                        ++fi;
                        if (fi == N) fi = 0;
                    }
                    if (last_unknown[fi] + S <= p) {
                        if (fkmer[fi] <= rkmer[fi]) {
                            mers.emplace_back(fkmer[fi],p);
                        } else {
                            mers.emplace_back(rkmer[fi],-1*p);
                        }
                    }
                }
            }
        }
        return false;
    }

    const bool next_element(std::string &rec, std::vector<MinPosIDX> &mers){
        for (auto ni=0; ni < N; ++ni) {
            fkmer[ni] =0;
            rkmer[ni] = 0;
            last_unknown[ni] = 0;
        }
        fi=0;
        uint64_t p(0);
        while (p < rec.size()) {
            //fkmer: grows from the right (LSB)
            //rkmer: grows from the left (MSB)
            bases++;
            for (auto ni=0;ni<N;++ni) {
                cycle_pos[ni]++;
                if (cycle_pos[ni]==N)cycle_pos[ni]=0;
                if (cycle_pos[ni]<M) {
                    fillKBuf(rec[p], fkmer[ni], rkmer[ni], last_unknown[ni]);
                }
                p++;
                //if we are at p, the skip-mer started at p-S is now done
                if (p>=S-1) {
                    if (p == S - 1) fi = 0;
                    else {
                        ++fi;
                        if (fi == N) fi = 0;
                    }
                    if (last_unknown[fi] + S <= p) {
                        if (fkmer[fi] <= rkmer[fi]) {
                            mers.emplace_back(fkmer[fi],p);
                        } else {
                            mers.emplace_back(rkmer[fi],-1*p);
                        }
                    }
                }
            }
        }
        return false;
    }

    // TODO: Adjust for when K is larger than what fits in uint64_t!
    const bool next_element(std::vector<MinPosIDX> &mers) {
        uint64_t p(0);
        while (p < currentRecord.seq.size()) {
            //fkmer: grows from the right (LSB)
            //rkmer: grows from the left (MSB)
            bases++;
            for (auto ni=0;ni<N;++ni) {
                cycle_pos[ni]++;
                if (cycle_pos[ni]==N)cycle_pos[ni]=0;
                if (cycle_pos[ni]<M) {
                    fillKBuf(currentRecord.seq[p], fkmer[ni], rkmer[ni], last_unknown[ni]);
                }
                p++;
                //if we are at p, the skip-mer started at p-S is now done
                if (p>=S-1) {
                    if (p == S - 1) fi = 0;
                    else {
                        ++fi;
                        if (fi == N) fi = 0;
                    }
                    if (last_unknown[fi] + S <= p) {
                        if (fkmer[fi] <= rkmer[fi]) {
                            mers.emplace_back(fkmer[fi],p);
                        } else {
                            mers.emplace_back(rkmer[fi],-1*p);
                        }
                    }
                }
            }
        }
        return false;
    }

private:
    FileRecord currentRecord;
    uint64_t bases;
};


#endif //BSG_SKIPMERINDEXFACTORY_HPP
