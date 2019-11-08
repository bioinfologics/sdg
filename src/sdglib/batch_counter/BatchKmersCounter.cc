//
// Created by Luis Yanes (EI) on 2019-06-28.
//

#include <atomic>
#include "BatchKmersCounter.hpp"
#include <sdglib/workspace/WorkSpace.hpp>

std::shared_ptr<KmerList>
BatchKmersCounter::kmerCountOMPDiskBased(uint8_t K, PairedReadsDatastore const &reads,
                                         unsigned minCount, std::string tmpdir, std::string workdir,
                                         unsigned char disk_batches) {

    //If size larger than batch (or still not enough cpus used, or whatever), Lauch 2 tasks to sort the 2 halves, with minFreq=0

    sdglib::OutputLog() << "disk-based kmer counting with "<<(int) disk_batches<<" batches"<<std::endl;
    uint64_t total_kmers_in_batches=0;
    for (auto batch=0;batch < disk_batches;batch++) {
        uint64_t to = (batch+1) * reads.size()/disk_batches;
        uint64_t from= batch * reads.size()/disk_batches;
        auto kmer_list = kmerCountOMP(K, reads,from,to);
        std::ofstream batch_file(tmpdir+"/kmer_count_batch_"+std::to_string((int)batch),std::ios::out | std::ios::trunc | std::ios::binary);
        batch_file.write((const char *)kmer_list->kmers,sizeof(KMerNodeFreq_s)*kmer_list->size);
        batch_file.close();
        sdglib::OutputLog() << "batch "<<(int) batch<<" done and dumped with "<<kmer_list->size<< " kmers" <<std::endl;
        total_kmers_in_batches+=kmer_list->size;
    }

    //now a multi-merge sort between all batch files into the Dict
    sdglib::OutputLog() << "merging from disk"<<std::endl;
    //open all batch files
    std::ifstream dbf[disk_batches];
    bool dbf_active[disk_batches];
    KMerNodeFreq_s next_knf_from_dbf[disk_batches];
    uint finished_files=0;
    for (auto i=0;i<disk_batches;i++){
        dbf[i].open(tmpdir+"/kmer_count_batch_"+std::to_string((int)i),std::ios::in | std::ios::binary);
        dbf[i].read((char *)&next_knf_from_dbf[i],sizeof(KMerNodeFreq_s));
        dbf_active[i]=true;
    }
    //set all finished flags to false
    //TODO: stupid minimum search
    KMerNodeFreq_s current_kmer;
    //while finished_files<batches
    bool first=true;
    uint min=0;
    for (auto i=1;i<disk_batches;++i)
        if (dbf_active[i]){
            if (next_knf_from_dbf[i]<next_knf_from_dbf[min]) min=i;
        }
    current_kmer=next_knf_from_dbf[min];
    current_kmer.count=0;
    uint64_t used = 0,not_used=0;
    uint64_t hist[256];
    for (auto &h:hist) h=0;
    std::shared_ptr<KmerList> kmerlist=std::make_shared<KmerList>();
    //size_t last_kmer=0;
    size_t alloc_block=10000000;
    while (finished_files<disk_batches) {
        //find minimum of the non-finished files
        uint min=0;
        while (!dbf_active[min])++min;
        for (auto i=1;i<disk_batches;++i)
            if (dbf_active[i]){
                if (next_knf_from_dbf[i]<next_knf_from_dbf[min]) min=i;
            }
        //larger than current kmer?
        if (next_knf_from_dbf[min] > current_kmer) {
            ++hist[std::min(255,(int)current_kmer.count)];
            if (current_kmer.count>=minCount) {
                //(*dict)->insertEntryNoLocking(BRQ_Entry((BRQ_Kmer) current_kmer, current_kmer.kc));
                //kmerlist.push_back(current_kmer);
                if (used>=kmerlist->size) kmerlist->resize(kmerlist->size+alloc_block);
                kmerlist->kmers[used]=current_kmer;
                ++used;
            }
            else ++not_used;
            current_kmer=next_knf_from_dbf[min];
        } else {
            current_kmer.combine(next_knf_from_dbf[min]);
        }
        //advance min file
        dbf[min].read((char *)&next_knf_from_dbf[min],sizeof(KMerNodeFreq_s));
        if ( dbf[min].eof() ) {
            dbf_active[min]=false;
            ++finished_files;
        }
    }

    ++hist[(int)current_kmer.count];
    if (current_kmer.count>=minCount) {
        //kmerlist.push_back(current_kmer);
        if (used>=kmerlist->size) kmerlist->resize(kmerlist->size+alloc_block);
        kmerlist->kmers[used]=current_kmer;
        used++;
    }
    else not_used++;
    kmerlist->resize(used);
    for (auto i=0;i<disk_batches;i++) {
        dbf[i].close();
        std::remove((tmpdir + "/kmer_count_batch_" +std::to_string((int)i)).c_str());
    }
    sdglib::OutputLog()<< used << "/" << used+not_used << " kmers with Freq >= " << minCount << std::endl;
    if (""!=workdir) {
        std::ofstream kff(workdir + "/small_K.freqs");
        for (auto i = 1; i < 256; i++) kff << i << ", " << hist[i] << std::endl;
        kff.close();
    }
    return kmerlist;
}

std::shared_ptr<KmerList>
BatchKmersCounter::kmerCountOMP(uint8_t K, PairedReadsDatastore const &reads, uint64_t gfrom, uint64_t gto, uint64_t batch_size) {
    uint64_t totalKmers(0);
    //Compute how many "batches" will be used. and malloc a structure for them and a bool array to mark them "ready to process".
    //optionally, there could be a total count of "ready kmers" to set a limit for memory usage, if total count is too high, slaves would wait
    if (batch_size == 0) batch_size = (gto - gfrom) / (4 * omp_get_max_threads()) + 1;
    const uint64_t batches = ((gto - gfrom) + batch_size - 1) / batch_size;
    sdglib::OutputLog() << "OMP-merge kmer counting in " << batches << " batches of " << batch_size << " reads" << std::endl;

    std::vector<std::atomic_uint_fast8_t *> batch_status; //level->batch->status

    std::vector<std::vector<std::shared_ptr<KmerList>>> batch_lists; //level->batch-> pointer to batch
    uint16_t levels = 0;
    for (auto elements = batches; elements > 1; elements = (elements + 1) / 2) {
        batch_status.push_back((std::atomic_uint_fast8_t *) calloc(sizeof(std::atomic_uint_fast8_t), elements));
        batch_lists.emplace_back();
        batch_lists.back().resize(elements);
//        sdglib::OutputLog() << "level " << levels << " created with " << elements << " elements " << std::endl;
        ++levels;
    }
    sdglib::OutputLog() << "Created " << levels << " levels" << std::endl;
    std::atomic_uint_fast8_t level_count[levels]; //level->count
    for (auto &l:level_count)l = 0;

#pragma omp parallel shared(reads) reduction(+:totalKmers)
    {
        ReadSequenceBuffer bprsg(reads,100000,1000);
        std::vector<__uint128_t> read_kmers;
#pragma omp for schedule(dynamic)
        for (auto batch = 0; batch < batches; ++batch) {
            //==== Part 1: actually counting ====
            uint64_t from = gfrom + batch * batch_size;
            uint64_t read_count = (batch < batches - 1) ? batch_size : (gto - gfrom) - batch_size * (batches - 1);
            uint64_t to = from + read_count;

            std::shared_ptr<KmerList> local_kmer_list = std::make_shared<KmerList>();
            uint64_t total_good_length = reads.readsize*(to-from+1);
            local_kmer_list->resize(total_good_length);
            //Populate the kmer list
            uint64_t last_kmer=0;

            // Replace with generating all the kmers from the reads!
            StreamKmerFactory128 skf(K);

            for (uint64_t rid=from+1; rid<= to; ++rid) {
                skf.produce_all_kmers(bprsg.get_read_sequence(rid), read_kmers);
                for (const auto &rk: read_kmers){
                    memcpy((char*) &local_kmer_list->kmers[last_kmer].kdata, (char *) &rk, 16);
                    local_kmer_list->kmers[last_kmer].count=1;
                    ++last_kmer;
                }
                read_kmers.clear();
            }

            totalKmers += last_kmer;
            //std::cout<<last_kmer<<" kmers inserted in a batch"<<std::endl;
            local_kmer_list->resize(last_kmer);
            local_kmer_list->sort();
            local_kmer_list->uniq();
            //std::cout<<local_kmer_list->size<<" unique kmers in a batch"<<std::endl;
            //std::sort(local_kmer_list->begin(), local_kmer_list->end());
            //collapse_entries(*local_kmer_list);

            //==== Part 2: merging till no results of same level available ====
            //merge_level=0
            uint16_t merge_level = 0;
            bool just_merged = true;
            //while (just merged)
            while (just_merged) {
                //   insert the list into this merge level's queue
                uint16_t slot = level_count[merge_level]++;
                if (slot == batch_lists[merge_level].size() - 1) {
#pragma omp critical
                    sdglib::OutputLog() << "level " << (merge_level+1) << " done." << std::endl;
                }
                //   if insertion number is odd or last batch:
                if (merge_level < levels - 1 and (slot % 2 == 1 or slot == batch_lists[merge_level].size() - 1)) {
                    // mix with the previous even number (using own list as base, this way we preserve locality)
                    if (slot % 2 == 1) {
                        while (batch_status[merge_level][slot - 1] != 1)
                            usleep(10); //wait for the previous batch to be finished.
                        batch_status[merge_level][slot - 1] = 2;
                        batch_status[merge_level][slot] = 2;
                        //inplace_count_merge(local_kmer_list, batch_lists[merge_level][slot - 1]);
                        local_kmer_list->merge(*batch_lists[merge_level][slot - 1]);
                        batch_lists[merge_level][slot - 1].reset();
                    }
                    //      increase level
                    ++merge_level;
                    just_merged = true;
                } else {
                    // no batch available for merging, next thread at this level will pick up
                    batch_lists[merge_level][slot] = local_kmer_list;
                    local_kmer_list.reset();
                    batch_status[merge_level][slot] = 1;
                    just_merged = false;
                }
                //IDEA: fixed partitioning can be used to keep 2*thread-number lists for longer.
            }
        }

    }


//    sdglib::OutputLog() << "Top level merge starting" << std::endl;
    //inplace_count_merge(batch_lists.back()[0], batch_lists.back()[1]);
    batch_lists.back()[0]->merge(*batch_lists.back()[1]);
//    sdglib::OutputLog() << "Top level merge done" << std::endl;
    std::shared_ptr<KmerList> kmer_list = batch_lists.back()[0];
    batch_lists.back()[0].reset();
    batch_lists.back()[1].reset();
    for (auto &bs:batch_status) free(bs);
//    sdglib::OutputLog() << "cleanup done" << std::endl;
    sdglib::OutputLog() << "Total kmers processed " << totalKmers << std::endl;
    return kmer_list;
}

std::shared_ptr<KmerList>
BatchKmersCounter::buildKMerCount(uint8_t K, PairedReadsDatastore const &reads, unsigned minCount,
                                  std::string workdir, std::string tmpdir, unsigned char disk_batches) {
    std::shared_ptr<KmerList> spectrum=std::make_shared<KmerList>();
    sdglib::OutputLog() << "creating kmers from reads..." << std::endl;
    if (1 >= disk_batches) {
        uint64_t hist[256]={0};
        uint64_t used=0;
        spectrum=kmerCountOMP(K, reads, 0, reads.size(), 0);
        auto witr=spectrum->kmers;
        auto send=spectrum->kmers+spectrum->size;
        for (auto itr=spectrum->kmers;itr!=send;++itr){
            ++hist[itr->count];
            if (itr->count>=minCount) {
                (*witr)=(*itr);
                ++witr;
                ++used;
            }

        }
        sdglib::OutputLog() << used << "/" << spectrum->size << " " << int(K)<<"-mers with Freq >= " << minCount << std::endl;
        spectrum->resize(used);
        std::ofstream kff(workdir + "/small_K.freqs");
        for (auto i = 1; i < 256; i++) kff << i << ", " << hist[i] << std::endl;
        kff.close();
        //todo: filter kmers to minCount
    } else {
        if ("" == tmpdir) tmpdir = workdir;
        spectrum=kmerCountOMPDiskBased(K, reads, minCount, tmpdir, workdir, disk_batches);
    }

    return spectrum;
}

void KmerList::sort() {
    std::sort(kmers,kmers+size);
}

void KmerList::uniq() {
    auto wptr = kmers;
    auto endptr = kmers + size;
    for (auto rptr = kmers; rptr < endptr; ++wptr) {
        *wptr = *rptr;
        ++rptr;
        while (rptr < endptr and *wptr == *rptr) {
            wptr->combine(*rptr);
            ++rptr;
        }
    }
    resize(wptr - kmers);
}

void KmerList::dump(std::string filename) {
    std::ofstream batch_file(filename, std::ios::out | std::ios::trunc | std::ios::binary);
    uint64_t total_kmers = size;
    batch_file.write((const char *) &total_kmers, sizeof(uint64_t));
    batch_file.write((const char *) kmers, sizeof(KMerNodeFreq_s) * size);
    batch_file.close();
}

void KmerList::load(std::string filename) {
    std::ifstream batch_file(filename, std::ios::in | std::ios::binary);
    uint64_t total_kmers;
    batch_file.read((char *) &total_kmers, sizeof(uint64_t));
    resize(total_kmers);
    sdglib::OutputLog() <<"Reading "<<size<<" kmers"<<std::endl;
    batch_file.read((char *) kmers, sizeof(KMerNodeFreq_s) * size);
    batch_file.close();
}

void KmerList::resize(size_t new_size) {
    if (size != new_size) {
//        std::cout << " allocating space for "<< new_size <<" elements: " << sizeof(KMerNodeFreq_s) * new_size <<std::endl;
        if (0==size) kmers = (KMerNodeFreq_s *) malloc(sizeof(KMerNodeFreq_s) * new_size);
        else kmers = (KMerNodeFreq_s *) realloc(kmers, sizeof(KMerNodeFreq_s) * new_size);
//        if (new_size>0 and kmers == nullptr) std::cout << " realloc error!!! "<<std::endl;
        size = new_size;
    }
}

KmerList::~KmerList() {
    clear();
}

void KmerList::clear() {
    if (nullptr != kmers and 0 != size) {
        free(kmers);
        size = 0;
        kmers = nullptr;
    }
}

void KmerList::merge(KmerList &other) {
    //std::cout <<"merging batches of size "<<size<<" and "<<other.size<<std::endl;
    auto itr1 = kmers, end1 = kmers + size;
    auto itr2 = other.kmers, end2 = other.kmers + other.size, itr2w = other.kmers;
    //std::cout<<"merging, accumulating on counts1"<<std::endl;
    while (itr2 != end2) {
        while (itr1 != end1 and *itr1 < *itr2) ++itr1;
        if (itr1 != end1 and *itr1 == *itr2) {
            //combine_Entries(*itr1,*itr2);
            itr1->combine(*itr2);
            ++itr1;
            ++itr2;
        }
        while (itr2 != end2 and (itr1 == end1 or *itr2 < *itr1)) {
            *itr2w = *itr2;
            ++itr2w;
            ++itr2;
        }
    }
    //shrink counts2 to the size of its remaining elements
    //std::cout<<" equal elements merge done, now resizing and merging "<< (itr2w - other.kmers) << " new elements "<<std::endl;
    //other.resize(itr2w - other.kmers);
    other.size=itr2w - other.kmers;
    //counts2->shrink_to_fit();
    //expand counts1 to allow the insertion of the unique values on counts2
    //std::cout<<"merging, resizing counts1 to " << counts1.size()+counts2.size() << std::endl;
    auto old_size=size;
    resize(size + other.size);
    auto ritr1 = kmers + old_size - 1;

    //merge-sort from the bottom into count1.
    //std::cout<<"merging, final merging"<<std::endl;
    auto writr1 = kmers + size - 1, rend1 = kmers;
    auto ritr2 = other.kmers + other.size - 1, rend2 = other.kmers;

    while (writr1 >= rend1) {
        if (ritr2 >= rend2 and (ritr1 < rend1 or *ritr2 > *ritr1)) {
            memcpy(writr1, ritr2, sizeof(KMerNodeFreq_s));
            --ritr2;
        } else {
            memcpy(writr1, ritr1, sizeof(KMerNodeFreq_s));
            --ritr1;
        }
        --writr1;
    }
    other.clear();
    //std::cout<<"merge done, new size "<<size<<std::endl;
}


std::vector<__uint128_t> BatchKmersCounter::countKmersToList(const PairedReadsDatastore &ds, int k, int min_coverage, int num_batches) {


    auto kmer_list = BatchKmersCounter::buildKMerCount(k, ds, min_coverage, ".", ".", num_batches);
    std::vector<__uint128_t> kmers;
    kmers.reserve(kmer_list->size);
    __uint128_t kmer;
    for (uint64_t i = 0; i < kmer_list->size; ++i) {
        memcpy(&kmer, &kmer_list->kmers[i], 16);
        kmers.emplace_back(kmer);
    }
    sdglib::sort(kmers.begin(),kmers.end());
    return kmers;
}