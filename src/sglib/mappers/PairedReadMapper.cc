//
// Created by Bernardo Clavijo (EI) on 03/11/2017.
//
#include <iostream>
#include <iomanip>
#include <cassert>
#include <atomic>
#include <omp.h>
#include "PairedReadMapper.hpp"
#include <parallel/algorithm>

void PairedReadMapper::write(std::ofstream &output_file) {
    //read-to-node
    uint64_t count=read_to_node.size();
    output_file.write((const char *) &count,sizeof(count));
    output_file.write((const char *) read_to_node.data(),sizeof(sgNodeID_t)*count);
    //mappings
    count=reads_in_node.size();
    output_file.write((const char *) &count,sizeof(count));
    for (auto i=0;i<count;++i) {
        uint64_t mcount=reads_in_node[i].size();
        output_file.write((const char *) &mcount,sizeof(mcount));
        output_file.write((const char *) reads_in_node[i].data(), sizeof(ReadMapper) * mcount);
    }
}

void PairedReadMapper::read(std::ifstream &input_file) {
    uint64_t count;
    input_file.read(( char *) &count,sizeof(count));
    read_to_node.resize(count);
    input_file.read(( char *) read_to_node.data(),sizeof(sgNodeID_t)*count);
    input_file.read(( char *) &count,sizeof(count));
    reads_in_node.resize(count);
    for (auto i=0;i<count;++i) {
        uint64_t mcount;
        input_file.read(( char *) &mcount,sizeof(mcount));
        reads_in_node[i].resize(mcount);
        input_file.read(( char *) reads_in_node[i].data(), sizeof(ReadMapper) * mcount);
    }
}

void PairedReadMapper::update_graph_index() {
    sglib::OutputLog(sglib::INFO) << "Indexing graph..."<<std::endl;
    const int k = 31;
    uint64_t total_k=0;
    std::vector<KmerIDX> kidxv;
    for (auto &n:sg.nodes) if (n.sequence.size()>=k) total_k+=n.sequence.size()+1-k;
    kidxv.reserve(total_k);
    FastaRecord r;
    kmerIDXFactory <FastaRecord> kcf({k});
    for (sgNodeID_t n=1;n<sg.nodes.size();++n){
        if (sg.nodes[n].sequence.size()>=k){
            r.id=n;
            r.seq=sg.nodes[n].sequence;
            kcf.setFileRecord(r);
            kcf.next_element(kidxv);
        }
    }
    sglib::OutputLog(sglib::INFO)<<kidxv.size()<<" kmers in total"<<std::endl;
    sglib::OutputLog(sglib::INFO) << "  Sorting..."<<std::endl;
    std::sort(kidxv.begin(),kidxv.end());
    sglib::OutputLog(sglib::INFO) << "  Merging..."<<std::endl;
    auto wi=kidxv.begin();
    auto ri=kidxv.begin();
    auto nri=kidxv.begin();
    while (ri<kidxv.end()){
        while (nri!=kidxv.end() and nri->kmer==ri->kmer) ++nri;
        if (nri-ri==1) {
            *wi=*ri;
            ++wi;
        }
        ri=nri;

    }

    kidxv.resize(wi-kidxv.begin());
    sglib::OutputLog(sglib::INFO)<<kidxv.size()<<" unique kmers in index, creating map"<<std::endl;
    std::unordered_set<int32_t> seen_contigs;
    seen_contigs.reserve(sg.nodes.size());
    for (auto &kidx :kidxv) {
        kmer_to_graphposition[kidx.kmer]={kidx.contigID,kidx.pos};
        seen_contigs.insert((kidx.contigID>0?kidx.contigID:-kidx.contigID));
    }
    sglib::OutputLog(sglib::INFO)<<seen_contigs.size()<<" nodes with indexed kmers"<<std::endl;
    reads_in_node.resize(sg.nodes.size());
}



class StreamKmerFactory : public  KMerFactory {
public:
    explicit StreamKmerFactory(uint8_t k) : KMerFactory(k){}
    inline void produce_all_kmers(const char * seq, std::vector<KmerIDX> &mers){
        // TODO: Adjust for when K is larger than what fits in uint64_t!
        last_unknown=0;
        fkmer=0;
        rkmer=0;
        auto s=seq;
        while (*s!='\0' and *s!='\n') {
            //fkmer: grows from the right (LSB)
            //rkmer: grows from the left (MSB)
            fillKBuf(*s, 0, fkmer, rkmer, last_unknown);
            if (last_unknown >= K) {
                if (fkmer <= rkmer) {
                    // Is fwd
                    mers.emplace_back(fkmer);
                    mers.back().contigID=1;
                } else {
                    // Is bwd
                    mers.emplace_back(rkmer);
                    mers.back().contigID=-1;
                }
            }
            ++s;
        }
    }
};

void PairedReadMapper::remap_all_reads() {
    update_graph_index();
    for (auto &rtn:read_to_node) rtn=0;
    for (auto &rin:reads_in_node) rin.clear();
    map_reads();
}

void PairedReadMapper::map_reads(const std::unordered_set<uint64_t> &reads_to_remap) {
    const int k = 31;
    std::atomic_int64_t nokmers(0);
    read_to_node.resize(datastore.size()+1);
    if (not reads_to_remap.empty())
        sglib::OutputLog()<<reads_to_remap.size()<<" selected reads / "<<read_to_node.size()-1<<" total"<<std::endl;

    /*
     * Read mapping in parallel,
     */
    uint64_t thread_mapped_count[omp_get_max_threads()],thread_total_count[omp_get_max_threads()],thread_multimap_count[omp_get_max_threads()];
    std::vector<ReadMapper> thread_mapping_results[omp_get_max_threads()];
    sglib::OutputLog(sglib::LogLevels::DEBUG)<<"Private mapping initialised for "<<omp_get_max_threads()<<" threads"<<std::endl;
#pragma omp parallel
    {
        const int min_matches=1;
        std::vector<KmerIDX> readkmers;
        StreamKmerFactory skf(31);
        ReadMapper mapping;
        auto blrs=BufferedPairedSequenceGetter(datastore,128*1024,260);
        auto & private_results=thread_mapping_results[omp_get_thread_num()];
        auto & mapped_count=thread_mapped_count[omp_get_thread_num()];
        auto & total_count=thread_total_count[omp_get_thread_num()];
        auto & multimap_count=thread_multimap_count[omp_get_thread_num()];
        mapped_count=0;
        total_count=0;
        multimap_count=0;
        bool c ;
        //std::cout<<omp_get_thread_num()<<std::endl;
#pragma omp for
        for (uint64_t readID=1;readID<read_to_node.size();++readID) {
            mapping.read_id = readID;
            //this enables partial read re-mapping by setting read_to_node to 0
            if ((reads_to_remap.size()>0 and reads_to_remap.count(mapping.read_id)>0) or (reads_to_remap.empty() and 0==read_to_node[mapping.read_id])) {
                mapping.node = 0;
                mapping.unique_matches = 0;
                mapping.first_pos = 0;
                mapping.last_pos = 0;
                mapping.rev = false;
                mapping.unique_matches = 0;
                //get all kmers from read
                auto seq=blrs.get_read_sequence(readID);
                readkmers.clear();
                skf.produce_all_kmers(seq,readkmers);
                if (readkmers.size()==0) {
                    ++nokmers;
                }
                for (auto &rk:readkmers) {
                    auto nk = kmer_to_graphposition.find(rk.kmer);
                    if (kmer_to_graphposition.end()!=nk) {
                        //get the node just as node
                        sgNodeID_t nknode = llabs(nk->second.node);
                        //TODO: sort out the sign/orientation representation
                        if (mapping.node == 0) {
                            mapping.node = nknode;
                            if ((nk->second.node > 0 and rk.contigID > 0) or
                                (nk->second.node < 0 and rk.contigID < 0))
                                mapping.rev = false;
                            else mapping.rev = true;
                            mapping.first_pos = nk->second.pos;
                            mapping.last_pos = nk->second.pos;
                            ++mapping.unique_matches;
                        } else {
                            //TODO:break mapping by change of direction and such
                            if (mapping.node != nknode) {
                                mapping.node = 0;
                                ++multimap_count;
                                break; //exit -> multi-mapping read! TODO: allow mapping to consecutive nodes
                            } else {
                                mapping.last_pos = nk->second.pos;
                                ++mapping.unique_matches;
                            }
                        }
                    }
                }
                if (mapping.node != 0 and mapping.unique_matches >= min_matches) {
                    //optimisation: just save the mapping in a thread private collection for now, have a single thread putting from that into de structure at the end
                    private_results.push_back(mapping);
                    ++mapped_count;
                }
            }
            auto tc = ++total_count;
            if (tc % 10000000 == 0) sglib::OutputLog()<< mapped_count << " / " << tc <<" ("<<multimap_count<<" multi-mapped)"<< std::endl;
        }
    }
    for (auto & tres:thread_mapping_results){
        //sglib::OutputLog(sglib::LogLevels::DEBUG)<<"mixing in "<<tres.size()<<" thread specific results"<<std::endl;
        for (auto &rm:tres){
            read_to_node[rm.read_id] = rm.node;
            reads_in_node[rm.node].emplace_back(rm);
        }
        tres.clear();
        tres.shrink_to_fit();
    }
    uint64_t mapped_count=0,total_count=0,multimap_count=0;
    for (auto i=0;i<omp_get_max_threads();++i){
        mapped_count+=thread_mapped_count[i];
        total_count+=thread_total_count[i];
        multimap_count+=thread_multimap_count[i];
    }
    sglib::OutputLog(sglib::LogLevels::INFO)<<"Reads without k-mers: "<<nokmers<<std::endl;
    sglib::OutputLog(sglib::LogLevels::INFO)<<"Reads mapped: "<<mapped_count<<" / "<<total_count<<" ("<<multimap_count<<" multi-mapped)"<<std::endl;
#pragma omp parallel for
    for (sgNodeID_t n=1;n<reads_in_node.size();++n){
        std::sort(reads_in_node[n].begin(),reads_in_node[n].end());
    }
}

///**
// * @brief Mapping of paired end read files.
// *
// * Reads are mapped through a unique k-mer index in process_reads_from_file.
// * R1 and R2 are processed independently. R1 gets odds ids, R2 gets the next even id, so [1,2] and [3,4] are pairs.
// *
// * @todo fix and test 10x tags processing
// * @todo enable some basic 10x tag statistics
// * @todo add support for LMP reads (i.e FR reads)
// * @todo add distribution computing on the fly
// * @todo support other kind of indexes and variable k-mer size
// */
//void PairedReadMapper::map_reads(std::string filename1, std::string filename2, prmReadType read_type, uint64_t max_mem) {
//    read1filename=std::move(filename1);
//    read2filename=std::move(filename2);
//    readType=read_type;
//    memlimit=max_mem;
//    remap_reads();
//}

void PairedReadMapper::remove_obsolete_mappings(){
    uint64_t nodes=0,reads=0;
    std::set<sgNodeID_t> updated_nodes;
    for (auto n=1;n<sg.nodes.size();++n) {
        if (sg.nodes[n].status==sgNodeDeleted) {
            updated_nodes.insert(n);
            updated_nodes.insert(-n);
            reads_in_node[n > 0 ? n : -n].clear();
            ++nodes;
        }
    }
    for (auto &read_node:read_to_node) {
        if (updated_nodes.count(read_node) !=0 ) {
            read_node=0;
            ++reads;
        }
    }
    std::cout << "obsolete mappings removed from "<<nodes<<" nodes, total "<<reads<<" reads."<<std::endl;
}

void PairedReadMapper::print_stats(){
    uint64_t none=0,single=0,both=0,same=0;
    for (uint64_t r1=1;r1<read_to_node.size();r1+=2){
        if (read_to_node[r1]==0) {
            if (read_to_node[r1+1]==0) ++none;
            else ++single;
        }
        else if (read_to_node[r1+1]==0) ++single;
        else {
            ++both;
            if (read_to_node[r1]==read_to_node[r1+1]) ++same;
        }
    }
    sglib::OutputLog()<<"Mapped pairs from "<<datastore.filename<<": None: "<<none<<"  Single: "<<single<<"  Both: "<<both<<" ("<<same<<" same)"<<std::endl;
}

std::vector<uint64_t> PairedReadMapper::size_distribution() {
    frdist.clear();
    frdist.resize(20000);
    rfdist.clear();
    rfdist.resize(20000);
    uint64_t frcount=0,rfcount=0;
    for (uint64_t r1=1;r1<read_to_node.size();r1+=2){

        if (read_to_node[r1]!=0 and read_to_node[r1]==read_to_node[r1+1]) {
            auto node=read_to_node[r1];
            ReadMapper rm1,rm2;
            for (auto rm:reads_in_node[node]){
                if (rm.read_id==r1) rm1=rm;
                if (rm.read_id==r1+1) rm2=rm;
            }
            if (rm1.first_pos>rm2.first_pos) std::swap(rm1,rm2);
            auto d=rm2.last_pos-rm1.first_pos;
            if (d>=rfdist.size()) continue;
            if (rm1.rev and !rm2.rev) {
                ++rfdist[d];
                ++rfcount;
            }
            if (!rm1.rev and rm2.rev) {
                ++frdist[d];
                ++frcount;
            }
        }
    }
    //std::cout<<"Read orientations:  FR: "<<frcount<<"  RF: "<<rfcount<<std::endl;
    if (frcount>rfcount){
        return frdist;
    } else return rfdist;
}