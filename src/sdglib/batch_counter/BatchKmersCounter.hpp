//
// Created by Luis Yanes (EI) on 2019-06-28.
//

#ifndef once

#include <string>
#include <memory>
#include <sdglib/datastores/PairedReadsDatastore.hpp>
#include <numeric>
#include <sdglib/datastores/ReadSequenceBuffer.hpp>

typedef struct __attribute__((__packed__)) KMerNodeFreq_s {
    uint64_t kdata[2];
    uint8_t count;
    inline const bool operator==(KMerNodeFreq_s const & other) const {
        //return 0==memcmp(&kdata,&other.kdata,2*sizeof(uint64_t));
        return (kdata[0]==other.kdata[0] and kdata[1]==other.kdata[1]);
    }
    inline const bool operator<(KMerNodeFreq_s const & other) const{
        //return -1==memcmp(&kdata,&other.kdata,2*sizeof(uint64_t));
        if (kdata[0]<other.kdata[0]) return true;
        if (kdata[0]>other.kdata[0]) return false;
        return kdata[1]<other.kdata[1];
    }
    inline const bool operator>(KMerNodeFreq_s const & other) const{
        if (kdata[0]>other.kdata[0]) return true;
        if (kdata[0]<other.kdata[0]) return false;
        return kdata[1]>other.kdata[1];
    }
    inline void combine(KMerNodeFreq_s const & other){
        uint16_t newcount=count+other.count;
        count = (newcount > UINT8_MAX) ? UINT8_MAX : newcount;
    }
    inline void merge(KMerNodeFreq_s const & other){
        uint16_t newcount=count+other.count;
        count = (newcount > UINT8_MAX) ? UINT8_MAX : newcount;
    }

};

class KmerList{
public:
    ~KmerList();
    void clear();
    void merge(KmerList & other);
    void sort();
    void uniq();
    void dump(std::string filename);
    void load(std::string filename);
    void resize(size_t new_size);
    KMerNodeFreq_s * kmers = nullptr;
    size_t size=0;
};

class BatchKmersCounter {

    static std::shared_ptr<KmerList> kmerCountOMP( uint8_t K, PairedReadsDatastore const& reads,
                                            uint64_t gfrom, uint64_t gto, uint64_t batch_size=0);

    static std::shared_ptr<KmerList> kmerCountOMP( uint8_t K, LongReadsDatastore const& reads,
                                            uint64_t gfrom, uint64_t gto, uint64_t batch_size=0);


    static std::shared_ptr<KmerList> kmerCountOMPDiskBased( uint8_t K, PairedReadsDatastore const& reads, unsigned minCount,
                                                     std::string tmpdir, std::string workdir, unsigned char disk_batches=8);

    static std::shared_ptr<KmerList> kmerCountOMPDiskBased( uint8_t K, LongReadsDatastore const& reads, unsigned minCount,
                                                     std::string tmpdir, std::string workdir, unsigned char disk_batches=8);

public:
    static std::shared_ptr<KmerList> buildKMerCount( uint8_t K, PairedReadsDatastore const& reads, unsigned minCount,
                                              std::string workdir, std::string tmpdir,
                                              unsigned char disk_batches );

    static std::shared_ptr<KmerList> buildKMerCount( uint8_t K, LongReadsDatastore const& reads, unsigned minCount,
                                                     std::string workdir, std::string tmpdir,
                                                     unsigned char disk_batches );

    static std::vector<__uint128_t> countKmersToList(const PairedReadsDatastore &ds, int k, int min_coverage, int num_batches);

    static std::vector<__uint128_t> countKmersToList(const LongReadsDatastore &ds, int k, int min_coverage, int num_batches);

};


#endif //SDG_BATCHKMERSCOUNTER_HPP
