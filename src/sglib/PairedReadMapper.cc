//
// Created by Bernardo Clavijo (EI) on 03/11/2017.
//

#include <iostream>
#include <iomanip>
#include <cassert>
#include <atomic>
#include "PairedReadMapper.hpp"


#define mem_limit 4
#define GB 1024*1024*1024


uint64_t PairedReadMapper::process_reads_from_file(uint8_t k, uint16_t min_matches, std::vector<KmerIDX> &unique_kmers, std::string filename, uint64_t offset , bool is_tagged=false) {
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
        ReadMapping mapping;
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
                if (is_tagged){
                    //std::cout << "Read name: " << read.name << std::endl;
                    std::string barcode = read.name.substr(read.name.size() - 16);
                    #pragma omp critical(add_mapped_tagged)
                    read_ids_to_tags[mapping.read_id] = barcode;
                }
                ++mapped_count;
            }
            auto tc=++total_count;
            if (tc % 100000 == 0) std::cout << mapped_count << " / " << tc << std::endl;
#pragma omp critical(read_record)
            c = fastqReader.next_record(read);
        }

    }
    // somehow for my test data with 700 reads, totak count is 834 for r2...
    std::cout<<"Reads mapped: "<<mapped_count<<" / "<<total_count<<std::endl;
    fastqReader.getSummaryStatistics();
    return total_count;
}

/*
 * Mapping of paired end reads.
 * Reads are mapped through a unique k-mer index. R1 and R2 are processed independently.
 * R1 gets odds ids, R2 gets the next even id, so [1,2] and [3,4] are pairs.
 */
void PairedReadMapper::map_reads(std::string filename1, std::string filename2, prmReadType read_type) {
    std::cout<<"Mapping "<<prmReadTypeDesc[read_type]<<" reads from "<<filename1<<" and "<<filename2<<std::endl;

        uint64_t maxmem(mem_limit * GB);
        std::cout << "maxmem=" << maxmem << std::endl;
        /*
         * Reads a fasta file and generates the set of kmers and which entry number on the fasta they belong to
         * along with their orientation in the form of a int32_t with + for fwd and - for rev, the output is filtered
         * by min < coverage <= max
         */

        const int k = 31;
        const int max_coverage = 1;
        const int min_matches = 1;
        const std::string output_prefix("smr_files");
        SMR<KmerIDX,
            kmerIDXFactory<FastaRecord>,
            GraphNodeReader<FastaRecord>,
            FastaRecord, GraphNodeReaderParams, KMerIDXFactoryParams> kmerIDX_SMR({1,sg}, {k}, maxmem, 0, max_coverage,
                                                                          output_prefix);

        std::vector<KmerIDX> unique_kmers;

        // Get the unique_kmers from the file
        std::cout << "Indexing assembly... " << std::endl;
        unique_kmers = kmerIDX_SMR.read_from_file(output_prefix);

        std::vector<uint64_t> uniqKmer_statistics(kmerIDX_SMR.summaryStatistics());
        std::cout << "Number of " << int(k) << "-kmers seen in assembly " << uniqKmer_statistics[0] << std::endl;
        //std::cout << "Number of " << int(k) << "-kmers that appear more than 0 < kmer_coverage <= " << max_coverage
        //          << " in assembly " << uniqKmer_statistics[1] << std::endl;
        std::cout << "Number of contigs from the assembly " << uniqKmer_statistics[2] << std::endl;
    if (read_type==prmPE) {
        auto r1c = process_reads_from_file(k, min_matches, unique_kmers, filename1, 1);
        auto r2c = process_reads_from_file(k, min_matches, unique_kmers, filename2, 2);
        //now populate the read_to_node array
        assert(r1c == r2c);
        read_to_node.resize(r1c * 2 + 1, 0);
        for (auto &rin:reads_in_node)
            for (auto &mr:rin)
                read_to_node[mr.read_id] = mr.node;

    } else if (read_type==prm10x){ // can no longer access read name here, but dont' know how big to make it otherwise
        auto r1c = process_reads_from_file(k, min_matches, unique_kmers, filename1, 1, true);
        auto r2c = process_reads_from_file(k, min_matches, unique_kmers, filename2, 2, true);
        //now populate the read_to_node array
        assert(r1c == r2c);
        read_to_node.resize(r1c * 2 + 1, 0);
        read_to_tag.resize(r1c * 2 + 1);
        for (auto &rin:reads_in_node)
            for (auto &mr:rin) {
                // this is fine,
                read_to_node[mr.read_id] = mr.node;
                //std::cout << "Node: " << mr.node <<std::endl;
                //std::cout << "read id: " << mr.read_id << std::endl;
                //std::cout << "Tag: " << read_ids_to_tags[mr.read_id]<< std::endl;
                //std::cout << "gfdgfd" << read_to_tag[mr.read_id]; // malloc error too
                //std::cout << "read_to_tag size: " << read_to_tag.size() << "read_to_node size: " << read_to_node.size() << std::endl;
                //read_to_tag[mr.read_id] = read_ids_to_tags[mr.read_id];
                //std::cout << "read to tag " << read_to_tag[mr.read_id];

            }

        }
}

void PairedReadMapper::print_stats() {
    enum pair_status {pair_status_same,pair_status_different,pair_status_onlyr1,pair_status_onlyr2,pair_status_unmapped};
    uint64_t status_counts[]={0,0,0,0,0};
    //Now some starts about read mapping
    for (uint64_t rid=1; rid<read_to_node.size(); rid+=2){
        if (0==read_to_node[rid]) {
            if (0==read_to_node[rid+1]) ++status_counts[pair_status_unmapped];
            else ++status_counts[pair_status_onlyr2];
        }
        else {
            if (0==read_to_node[rid+1]) ++status_counts[pair_status_onlyr1];
            else {
                if (read_to_node[rid+1]==read_to_node[rid]) ++status_counts[pair_status_same];
                else ++status_counts[pair_status_different];
            }

        }
    }
    std::cout<<"---Pair mapping Stats ---"<<std::endl;
    std::cout<<"Both unmapped:   "<<status_counts[pair_status_unmapped]<<std::endl;
    std::cout<<"R1 only:         "<<status_counts[pair_status_onlyr1]<<std::endl;
    std::cout<<"R2 only:         "<<status_counts[pair_status_onlyr2]<<std::endl;
    std::cout<<"Both different:  "<<status_counts[pair_status_different]<<std::endl;
    std::cout<<"Both same :      "<<status_counts[pair_status_same]<<std::endl;
    std::cout<<"TOTAL     :      "<<read_to_node.size()<<std::endl;
    std::cout<<"---Node occupancy histogram ---"<<std::endl;
    uint64_t readcount[12];
    for (auto &rc:readcount)rc=0;
    for (sgNodeID_t n; n<reads_in_node.size();++n){
        auto c=reads_in_node[n].size();
        if (c>0) c=c/100+1;
        if (c>11) c=11;
        ++readcount[c];
    }
    std::cout<<"0: "<<readcount[0]<<std::endl;
    for (int i=1;i<11;++i) if (readcount[i]) std::cout<<(i>1? (i-1)*100: 1)<<"-"<<i*100-1<<": "<<readcount[i]<<std::endl;
    std::cout<<"1000+: "<<readcount[11]<<std::endl;

}

void PairedReadMapper::save_to_disk(std::string filename) {
    std::ofstream of(filename);
    //read-to-node
    uint64_t count=read_to_node.size();
    of.write((const char *) &count,sizeof(count));
    of.write((const char *) read_to_node.data(),sizeof(sgNodeID_t)*count);
    //reads-in-node
    count=reads_in_node.size();
    of.write((const char *) &count, sizeof(count));
    for (auto &rtn:reads_in_node) {
        count=rtn.size();
        of.write((const char *) &count, sizeof(count));
        of.write((const char *) rtn.data(),sizeof(ReadMapping)*count);
    }
}

void PairedReadMapper::load_from_disk(std::string filename) {

    std::ifstream inf(filename);
    uint64_t count;

    //read-to-node

    inf.read((char *) &count,sizeof(count));
    read_to_node.resize(count);
    inf.read((char *) read_to_node.data(),sizeof(sgNodeID_t)*count);
    //reads-in-node
    inf.read((char *) &count,sizeof(count));
    reads_in_node.resize(count);
    for (auto &rtn:reads_in_node) {
        inf.read((char *) &count, sizeof(count));
        rtn.resize(count);
        inf.read((char *) rtn.data(),sizeof(ReadMapping)*count);
    }

}

