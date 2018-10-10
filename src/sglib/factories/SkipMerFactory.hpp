//
// Created by Luis Yanes (EI) on 20/10/2017.
//

#ifndef SEQ_SORTER_SKIPKMERFACTORY_H
#define SEQ_SORTER_SKIPKMERFACTORY_H

#include <cstdint>
#include <iostream>
#include <vector>
#include <sglib/readers/FileReader.hpp>
#include <sglib/types/KmerTypes.hpp>

#define unlikely(x)     __builtin_expect((x),0)


class SkipMerFactory {
public:
    const uint8_t K;
    const uint8_t M;
    const uint8_t N;
    const int S;
    std::vector<int64_t> last_unknown{};
    std::vector<uint64_t> fkmer{};
    std::vector<uint64_t> rkmer{};
    std::vector<uint8_t> cycle_pos{};
    uint8_t fi{0};

    inline void fillKBuf(const char b, uint64_t &rkmer, uint64_t &fkmer, int64_t &last) {
        if (unlikely(b2f[b] == 4)) {
            last = 0;
            fkmer = ((fkmer << 2) + 0) & KMER_MASK;
            rkmer = (rkmer >> 2) + (((uint64_t) 3) << KMER_FIRSTOFFSET);
        } else {
            fkmer = ((fkmer << 2) + b2f[b]) & KMER_MASK;
            rkmer = (rkmer >> 2) + (((uint64_t) b2r[b]) << KMER_FIRSTOFFSET);
            ++last;
        }
    }

protected:
    explicit SkipMerFactory(uint8_t k, uint8_t m, uint8_t n) : K(k), M(m), N(n), KMER_FIRSTOFFSET((uint64_t) (k - 1) * 2),
                                      KMER_MASK((((uint64_t) 1) << (k * 2)) - 1), S(k+((k-1)/(int)m)*((int)n-(int)m)) {
        for (int i = 0; i < 255; i++) {
            b2f[i] = 4;
            b2r[i] = 4;
        }
        b2f['a'] = b2f['A'] = 0;
        b2f['c'] = b2f['C'] = 1;
        b2f['g'] = b2f['G'] = 2;
        b2f['t'] = b2f['T'] = 3;

        b2r['a'] = b2r['A'] = 3;
        b2r['c'] = b2r['C'] = 2;
        b2r['g'] = b2r['G'] = 1;
        b2r['t'] = b2r['T'] = 0;

        fkmer = std::vector<uint64_t>(n);
        rkmer = std::vector<uint64_t>(n);
        last_unknown = std::vector<int64_t>(n);
        cycle_pos = std::vector<uint8_t>(n);
    }

private:
    const uint64_t KMER_MASK;
    const uint64_t KMER_FIRSTOFFSET;
    char b2f[255]{4};
    char b2r[255]{4};
};

class SkipmerIndexFactory : protected SkipMerFactory {
public:
    SkipmerIndexFactory(uint8_t k, uint8_t m, uint8_t n) : SkipMerFactory(k, m, n) {}

    void setFileRecord(FastaRecord &rec) {
        currentRecord = rec;
        for (auto ni=0; ni < N; ++ni) {
            fkmer[ni] =0;
            rkmer[ni] = 0;
            last_unknown[ni] = 0;
        }
        fi=0;
    }
    // TODO: Adjust for when K is larger than what fits in uint64_t!
    const bool next_element(std::vector<std::pair<uint64_t,graphPosition>> &mers){
        for (auto ni=0; ni < N; ++ni) {
            fkmer[ni] =0;
            rkmer[ni] = 0;
            last_unknown[ni] = 0;
        }
        fi=0;
        uint64_t p(0);
        graphPosition pos;
        while (p < currentRecord.seq.size()) {
            //fkmer: grows from the right (LSB)
            //rkmer: grows from the left (MSB)
            for (auto ni=0;ni<N;++ni) {
                cycle_pos[ni]++;
                if (cycle_pos[ni]==N)cycle_pos[ni]=0;
                if (cycle_pos[ni]<M) {
                    fillKBuf(currentRecord.seq[p], fkmer[ni], rkmer[ni], last_unknown[ni]);
                }
                p++;
                //if we are at p, the skip-mer started at p-S is now done

                if (p < S-1) continue;

                if (p == S - 1) fi = 0;
                else {
                    ++fi;
                    if (fi == N) fi = 0;
                }
                if (unlikely(last_unknown[fi] + S > p)) continue;

                if (fkmer[fi] <= rkmer[fi]) {
                    pos.node=currentRecord.id;
                    pos.pos=p;
                    mers.emplace_back(fkmer[fi],pos);
                } else {
                    pos.node=-currentRecord.id;
                    pos.pos=p;
                    mers.emplace_back(rkmer[fi],pos);
                }

            }
        }
        return false;
    }



private:
    FastaRecord currentRecord;
    uint64_t bases;
};


#endif //SEQ_SORTER_SKIPKMERFACTORY_H
