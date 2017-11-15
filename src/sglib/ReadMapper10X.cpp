//
// Created by Katie Barr (EI) on 15/11/2017.
//

#include "ReadMapper10X.h"

uint64_t ReadMapper10X::process_reads_from_file(uint8_t k, uint16_t min_matches, std::vector<KmerIDX> &unique_kmers, std::string filename, uint64_t offset ) {
    std::cout<<"mapping reads!!!"<<std::endl;
    /*
     * Read mapping in parallel,
     */
    FastqReader<FastqRecord> fastqReader({0},filename);
    std::atomic<uint64_t> mapped_count(0),total_count(0);
#pragma omp parallel shared(fastqReader,reads_in_node)
    {
        FastqRecord read;
        std::vector<KmerIDX> readkmers;
        kmerIDXFactory<FastqRecord> kf({k});
        ReadMapping10X mapping;
        bool c ;
#pragma omp critical(read_record)
        c = fastqReader.next_record(read);
        while (c) {

            //get all kmers from read
            readkmers.clear();
            kf.setFileRecord(read);
            kf.next_element(readkmers);

            mapping.node = 0;
            mapping.unique_matches = 0;
            for (auto &rk:readkmers) {
                auto nk = std::lower_bound(unique_kmers.begin(), unique_kmers.end(), rk);
                if (nk->kmer == rk.kmer) {
                    //get the node just as node
                    sgNodeID_t nknode = (nk->contigID > 0 ? nk->contigID : -nk->contigID);
                    //TODO: sort out the sign/orientation representation
                    if (mapping.node == 0) {
                        mapping.node = nknode;
                        if ((nk->contigID > 0 and rk.contigID > 0) or (nk->contigID < 0 and rk.contigID < 0)) mapping.rev=false;
                        else mapping.rev=true;
                        mapping.first_pos = nk->pos;
                        mapping.last_pos = nk->pos;
                        ++mapping.unique_matches;
                    } else {
                        if (mapping.node != nknode) {
                            mapping.node = 0;
                            break; //exit -> multi-mapping read! TODO: allow mapping to consecutive nodes
                        } else {
                            mapping.last_pos = nk->pos;
                            ++mapping.unique_matches;
                        }
                    }
                }
            }
            if (mapping.node != 0 and mapping.unique_matches >= min_matches) {
                //TODO: set read id and add to map collection
                mapping.read_id=(read.id)*2+offset;
#pragma omp critical(add_mapped)
                reads_in_node[mapping.node].emplace_back(mapping);
                std::string barcode = read.name.substr(read.name.size() - 16);
                mapping.barcode = barcode;

                ++mapped_count;
            }
            auto tc=++total_count;
            if (tc % 100000 == 0) std::cout << mapped_count << " / " << tc << std::endl;
#pragma omp critical(read_record)
            c = fastqReader.next_record(read);
        }

    }
    std::cout<<"Reads mapped: "<<mapped_count<<" / "<<total_count<<std::endl;
    return total_count;
}