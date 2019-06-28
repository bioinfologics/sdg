//
// Created by Bernardo Clavijo (EI) on 2019-06-25.
//

#include "KmerCountsDatastore.hpp"

template<class T>
void add_count_to_kds( KmerCountsDatastore & kds, const std::string & count_name, const T & datastore){
    kds.count_names.emplace_back(count_name);
    kds.counts.emplace_back(kds.kindex.size());
    uint64_t present(0), absent(0), rp(0);
    sdglib::OutputLog(sdglib::INFO)<<"Populating lookup map"<<std::endl;
    auto kmer_map=kds.get_kindex_map();
    sdglib::OutputLog(sdglib::INFO)<<"Map populated, counting from datastore: " << datastore.filename << std::endl;

#pragma omp parallel
    {
        ReadSequenceBuffer bpsg(datastore);
        uint64_t thread_present(0), thread_absent(0), thread_rp(0);
        const size_t local_kmers_size = 2000000;
        std::vector<uint64_t> found_kmers;
        found_kmers.reserve(local_kmers_size);
        std::vector<KmerCount> readkmers;
        CStringKMerFactory cskf(31);
#pragma omp for schedule(static,10000)
        for (uint64_t rid = 1; rid <= datastore.size(); ++rid) {
            readkmers.clear();
            cskf.create_kmercounts(readkmers, bpsg.get_read_sequence(rid));

            for (auto &rk:readkmers) {
                auto findk = kmer_map.find(rk.kmer);
                if (kmer_map.end() != findk) {
                    //++thread_counts[findk->second];
                    found_kmers.emplace_back(findk->second);
                    if (found_kmers.size() == local_kmers_size) {
                        bool printstatus = false;
#pragma omp critical(results_merge)
                        {
                            auto &arc=kds.counts.back();
                            for (auto &x:found_kmers) if (arc[x] < UINT16_MAX) ++arc[x];
                            if (rp / 1000000 != (rp + thread_rp) / 1000000) printstatus = true;
                            rp += thread_rp;
                            present += thread_present;
                            absent += thread_absent;
                        }
                        found_kmers.clear();
                        thread_absent = 0;
                        thread_present = 0;
                        thread_rp = 0;
                        if (printstatus)
                            sdglib::OutputLog(sdglib::INFO) << rp << " reads processed " << present << " / "
                                                            << present + absent << " kmers found" << std::endl;
                    }
                    ++thread_present;
                } else ++thread_absent;


            }
            ++thread_rp;
        }
    }
    sdglib::OutputLog(sdglib::INFO) << rp << " reads processed " << present << " / " << present + absent
                                    << " kmers found" << std::endl;
}

void KmerCountsDatastore::index_sdg() {
    //add all k-mers from SDG
    counts.clear();
    count_names.clear();
    uint64_t t=0;
    for(auto &n:ws.sdg.nodes) if (n.sequence.size()>=k) t+=n.sequence.size()+1-k;
    kindex.reserve(t);
    StringKMerFactory skf(k);
    for(auto &n:ws.sdg.nodes) if (n.sequence.size()>=k) skf.create_kmers(n.sequence,kindex);
    //sort
    std::sort(kindex.begin(),kindex.end());

    //create a first count
    counts.emplace_back();
    count_names.emplace_back("sdg");
    auto &c=counts.back();
    c.reserve(kindex.size());

    //collapse, but save coverage to the first count
    auto wi=kindex.begin();
    auto ri=kindex.begin();
    for (;ri<kindex.end();++wi){
        *wi=*ri;
        c.emplace_back(1);
        while(++ri<kindex.end() and *ri==*wi) ++(c.back());
    }
    kindex.resize(c.size());
}

std::unordered_map<uint64_t,uint64_t> & KmerCountsDatastore::get_kindex_map() const {
    std::unordered_map<uint64_t,uint64_t> kmer_map;
    kmer_map.reserve(kindex.size());
    for (uint64_t i=0;i<kindex.size();++i) kmer_map[kindex[i]]=i;
    return kmer_map;
}

void KmerCountsDatastore::add_count(const std::string & count_name, const PairedReadsDatastore & datastore){
    add_count_to_kds(*this,count_name,datastore);
}
void KmerCountsDatastore::add_count(const std::string & count_name, const LinkedReadsDatastore & datastore){
    add_count_to_kds(*this,count_name,datastore);
}
void KmerCountsDatastore::add_count(const std::string & count_name, const LongReadsDatastore & datastore){
    add_count_to_kds(*this,count_name,datastore);
}