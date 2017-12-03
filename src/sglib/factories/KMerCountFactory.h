//
// Created by Bernardo J. Clavijo (EI) on 3/12/2017.
//

#ifndef SMR_KMERCOUNTFACTORY_H
#define SMR_KMERCOUNTFACTORY_H

#include <vector>
#include <limits>
#include <tuple>
#include "sglib/factories/KMerFactory.h"
struct KmerCountFactoryParams {
    uint8_t k;
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
        if (other.count < 255-count) {
            count += other.count;
        }
        else count=255;
    }

    KmerCount max() {
        return {};
    }

    friend std::ostream& operator<<(std::ostream& os, const KmerCount& kmer) {
        os << kmer.kmer << "\t" << kmer.count;
        return os;
    }

    friend std::istream& operator>>(std::istream& is, const KmerCount& kmer) {
        is.read((char*)&kmer, sizeof(kmer));
        return is;
    }
    uint64_t kmer;
    uint8_t count;
};

namespace std {
    template <>
    struct hash<KmerCount> {
        size_t operator()(const KmerCount& k) const {
            return k.kmer;
        }
    };
}

template<typename FileRecord>
class KmerCountFactory : protected KMerFactory {
public:
    explicit KmerCountFactory(KmerCountFactoryParams params) : KMerFactory(params.k) {}

    ~KmerCountFactory() {
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
    const bool next_element(std::vector<KmerCount> &mers) {
        uint64_t p(0);
        while (p < currentRecord.seq.size()) {
            //fkmer: grows from the right (LSB)
            //rkmer: grows from the left (MSB)
            bases++;
            fillKBuf(currentRecord.seq[p], p, fkmer, rkmer, last_unknown);
            p++;
            if (last_unknown >= K) {
                if (fkmer <= rkmer) {
                    // Is fwd
                    mers.emplace_back(fkmer, 1);
                } else {
                    // Is bwd
                    mers.emplace_back(rkmer, 1);
                }
            }
        }
        return false;
    }

private:
    FileRecord currentRecord;
    uint64_t bases;
};
#endif //SMR_KMERCOUNTFACTORY_H