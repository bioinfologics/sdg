//
// Created by Bernardo Clavijo (EI) on 03/11/2017.
//
#include <iomanip>
#include <cassert>
#include <atomic>
#include <sdglib/utilities/omp_safe.hpp>
#include <fstream>
#include <iostream>
#include "PairedReadsMapper.hpp"
#include <sdglib/workspace/WorkSpace.hpp>
#include <sdglib/utilities/io_helpers.hpp>

const bsgVersion_t PairedReadsMapper::min_compat = 0x0001;

PairedReadsMapper::PairedReadsMapper(const WorkSpace &_ws, PairedReadsDatastore &_datastore) :
        ws(_ws),
        datastore(_datastore)
{
    reads_in_node.resize(ws.sdg.nodes.size());
}

void PairedReadsMapper::write(std::ofstream &output_file) {
    //read-to-node
    output_file.write((char *) &BSG_MAGIC, sizeof(BSG_MAGIC));
    output_file.write((char *) &BSG_VN, sizeof(BSG_VN));
    BSG_FILETYPE type(PairedMap_FT);
    output_file.write((char *) &type, sizeof(type));

    sdglib::write_flat_vector(output_file, read_to_node);
    //mappings
    sdglib::write_flat_vectorvector(output_file, reads_in_node);
}

void PairedReadsMapper::read(std::ifstream &input_file) {
    bsgMagic_t magic;
    bsgVersion_t version;
    BSG_FILETYPE type;
    input_file.read((char *) &magic, sizeof(magic));
    input_file.read((char *) &version, sizeof(version));
    input_file.read((char *) &type, sizeof(type));

    if (magic != BSG_MAGIC) {
        throw std::runtime_error("PairedReadsMapper file appears to be corrupted");
    }

    if (version < min_compat) {
        throw std::runtime_error("PairedReadsMapper file Incompatible version");
    }

    if (type != PairedMap_FT) {
        throw std::runtime_error("PairedReadsMapper file Incompatible file type");
    }

    sdglib::read_flat_vector(input_file, read_to_node);

    sdglib::read_flat_vectorvector(input_file, reads_in_node);

    populate_orientation();
}

void PairedReadsMapper::remap_all_reads() {
    for (auto &rtn:read_to_node) rtn=0;
    for (auto &rin:reads_in_node) rin.clear();
    map_reads();
}

void PairedReadsMapper::map_reads(const std::unordered_set<uint64_t> &reads_to_remap) {
    const int k = 31;
    std::atomic<int64_t> nokmers(0);
    reads_in_node.resize(ws.sdg.nodes.size());
    read_to_node.resize(datastore.size()+1);
    if (not reads_to_remap.empty())
        sdglib::OutputLog()<<reads_to_remap.size()<<" selected reads / "<<read_to_node.size()-1<<" total"<<std::endl;

    /*
     * Read mapping in parallel,
     */
    uint64_t thread_mapped_count[omp_get_max_threads()],thread_total_count[omp_get_max_threads()],thread_multimap_count[omp_get_max_threads()];
    std::vector<ReadMapping> thread_mapping_results[omp_get_max_threads()];
    sdglib::OutputLog(sdglib::LogLevels::DEBUG)<<"Private mapping initialised for "<<omp_get_max_threads()<<" threads"<<std::endl;
#pragma omp parallel
    {
        const int min_matches=1;
        std::vector<KmerIDX> readkmers;
        StreamKmerIDXFactory skf(31);
        ReadMapping mapping;
        auto blrs=BufferedPairedSequenceGetter(datastore,128*1024,datastore.readsize*2+2);
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
                    auto nk = ws.uniqueKmerIndex.find(rk.kmer);
                    if (ws.uniqueKmerIndex.end()!=nk) {
                        //get the node just as node
                        sgNodeID_t nknode = llabs(nk->second.node);
                        //TODO: sort out the sign/orientation representation
                        if (mapping.node == 0) {
                            mapping.node = nknode;
                            if ((nk->second.node > 0 and rk.contigID > 0) or
                                (nk->second.node < 0 and rk.contigID < 0))
                                mapping.rev = false;
                            else
                                mapping.rev = true;
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
            if (tc % 10000000 == 0) sdglib::OutputLog()<< mapped_count << " / " << tc <<" ("<<multimap_count<<" multi-mapped)"<< std::endl;
        }
    }
    for (auto & tres:thread_mapping_results){
        //sdglib::OutputLog(sdglib::LogLevels::DEBUG)<<"mixing in "<<tres.size()<<" thread specific results"<<std::endl;
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
    sdglib::OutputLog(sdglib::LogLevels::INFO)<<"Reads without k-mers: "<<nokmers<<std::endl;
    sdglib::OutputLog(sdglib::LogLevels::INFO)<<"Reads mapped: "<<mapped_count<<" / "<<total_count<<" ("<<multimap_count<<" multi-mapped)"<<std::endl;
#pragma omp parallel for
    for (sgNodeID_t n=1;n<reads_in_node.size();++n){
        std::sort(reads_in_node[n].begin(),reads_in_node[n].end());
    }
    populate_orientation();
}

void PairedReadsMapper::remap_all_reads63() {
    for (auto &rtn:read_to_node) rtn=0;
    for (auto &rin:reads_in_node) rin.clear();
    map_reads63();
}

void PairedReadsMapper::map_reads63(const std::unordered_set<uint64_t> &reads_to_remap) {
    const int k = 63;
    std::atomic<int64_t> nokmers(0);
    reads_in_node.resize(ws.sdg.nodes.size());
    read_to_node.resize(datastore.size()+1);
    if (not reads_to_remap.empty())
        sdglib::OutputLog()<<reads_to_remap.size()<<" selected reads / "<<read_to_node.size()-1<<" total"<<std::endl;

    /*
     * Read mapping in parallel,
     */
    uint64_t thread_mapped_count[omp_get_max_threads()],thread_total_count[omp_get_max_threads()],thread_multimap_count[omp_get_max_threads()];
    std::vector<ReadMapping> thread_mapping_results[omp_get_max_threads()];
    sdglib::OutputLog(sdglib::LogLevels::DEBUG)<<"Private mapping initialised for "<<omp_get_max_threads()<<" threads"<<std::endl;
#pragma omp parallel
    {
        const int min_matches=1;
        std::vector<KmerIDX128> readkmers;
        StreamKmerIDXFactory128 skf(63);
        ReadMapping mapping;
        auto blrs=BufferedPairedSequenceGetter(datastore,128*1024,datastore.readsize*2+2);
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
                    auto nk = ws.unique63merIndex.find(rk.kmer);
                    if (ws.unique63merIndex.end()!=nk) {
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
            if (tc % 10000000 == 0) sdglib::OutputLog()<< mapped_count << " / " << tc <<" ("<<multimap_count<<" multi-mapped)"<< std::endl;
        }
    }
    for (auto & tres:thread_mapping_results){
        //sdglib::OutputLog(sdglib::LogLevels::DEBUG)<<"mixing in "<<tres.size()<<" thread specific results"<<std::endl;
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
    sdglib::OutputLog(sdglib::LogLevels::INFO)<<"Reads without k-mers: "<<nokmers<<std::endl;
    sdglib::OutputLog(sdglib::LogLevels::INFO)<<"Reads mapped: "<<mapped_count<<" / "<<total_count<<" ("<<multimap_count<<" multi-mapped)"<<std::endl;
#pragma omp parallel for
    for (sgNodeID_t n=1;n<reads_in_node.size();++n){
        std::sort(reads_in_node[n].begin(),reads_in_node[n].end());
    }
    populate_orientation();
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

void PairedReadsMapper::remove_obsolete_mappings(){
    uint64_t nodes=0,reads=0;
    std::set<sgNodeID_t> updated_nodes;
    for (auto n=1;n<ws.sdg.nodes.size();++n) {
        if (ws.sdg.nodes[n].status==sgNodeDeleted) {
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
    sdglib::OutputLog(sdglib::INFO, false) << "obsolete mappings removed from "<<nodes<<" nodes, total "<<reads<<" reads."<<std::endl;
}

void PairedReadsMapper::print_status(){
    uint64_t none=0,single=0,both=0,same=0;
    for (uint64_t r1=1;r1<read_to_node.size();r1+=2){
        if (read_to_node[r1]==0) {
            if (read_to_node[r1+1]==0) ++none;
            else ++single;
        }
        else if (read_to_node[r1+1]==0) ++single;
        else {
            ++both;
            if (abs(read_to_node[r1])==abs(read_to_node[r1+1])) ++same;
        }
    }
    sdglib::OutputLog()<<"Mapped pairs from "<<datastore.filename<<": None: "<<none<<"  Single: "<<single<<"  Both: "<<both<<" ("<<same<<" same)"<<std::endl;
}

std::vector<uint64_t> PairedReadsMapper::size_distribution() {
    frdist.clear();
    frdist.resize(70000);
    rfdist.clear();
    rfdist.resize(70000);
    uint64_t frcount=0,rfcount=0;
    std::vector<int32_t> read_firstpos(read_to_node.size()),read_lastpos(read_to_node.size());
    std::vector<bool> read_rev(read_to_node.size());
    for (auto n=1;n<ws.sdg.nodes.size();++n) {
        for (auto &rm:reads_in_node[n]) {
            read_firstpos[rm.read_id]=rm.first_pos;
            read_lastpos[rm.read_id]=rm.last_pos;
            read_rev[rm.read_id]=rm.rev;
        }
    }
    for (uint64_t r1=1;r1<read_to_node.size();r1+=2){
        if (read_to_node[r1]!=0 and read_to_node[r1]==read_to_node[r1+1]) {
            auto node=read_to_node[r1];
            ReadMapping rm1,rm2;
            rm1.first_pos=read_firstpos[r1];
            rm1.last_pos=read_lastpos[r1];
            rm1.rev=read_rev[r1];
            rm2.first_pos=read_firstpos[r1+1];
            rm2.last_pos=read_lastpos[r1+1];
            rm2.rev=read_rev[r1+1];
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
    std::cout<<"Read orientations:  FR: "<<frcount<<"  RF: "<<rfcount<<std::endl;
    if (frcount>rfcount){
        return frdist;
    } else return rfdist;
}

void PairedReadsMapper::populate_orientation() {
    read_direction_in_node.clear();
    read_direction_in_node.resize(read_to_node.size());
    for (auto & nreads:reads_in_node){
        for (auto &rm:nreads){
            read_direction_in_node[rm.read_id]=rm.rev;
        }
    }
}

PairedReadsMapper& PairedReadsMapper::operator=(const PairedReadsMapper &other) {
    if (&ws.sdg != &other.ws.sdg and &datastore != &other.datastore) { throw std::runtime_error("Can only copy PairedReadsMapper from the same SequenceDistanceGraph and Datastore"); }
    if (&other == this) {
        return *this;
    }
    reads_in_node = other.reads_in_node;
    read_to_node = other.read_to_node;
    return *this;
}

PairedReadConnectivityDetail::PairedReadConnectivityDetail(const PairedReadsMapper &prm, sgNodeID_t source,
                                                           sgNodeID_t dest) {
    /*std::set<uint64_t> connecting_reads_s;
    std::set<uint64_t> connecting_reads_d;
    for (auto rm:prm.reads_in_node[llabs(source)]){
        auto rs=rm.read_id;
        auto rd=rs;
        if (rs%2==1) rd=rs+1;
        else rd=rs-1;
        connecting_reads_s.insert(rs);
        connecting_reads_d.insert(rd);
    }*/
    sgNodeID_t us=llabs(source);
    sgNodeID_t  ud=llabs(dest);
    for (auto rm:prm.reads_in_node[us]){
        uint64_t r1,r2;
        if (rm.read_id%2==1){
            r1=rm.read_id;
            r2=r1+1;
        }
        else{
            r1=rm.read_id-1;
            r2=r1+1;
        }
        if (prm.read_to_node[r1]==us){
            if (prm.read_to_node[r2]==ud){
                ++orientation_paircount[(prm.read_direction_in_node[r1]? 0:1)+(prm.read_direction_in_node[r2]? 0:2)];
            }
        }
        else if (prm.read_to_node[r1]==ud) {
            if (prm.read_to_node[r2]==us) {
                ++orientation_paircount[(prm.read_direction_in_node[r2] ? 0 : 1)+(prm.read_direction_in_node[r1] ? 0 : 2)];
            }

        }
    }

}

std::vector<uint64_t> PairedReadsMapper::get_node_readpairs_ids(sgNodeID_t nodeID) {
    std::vector<uint64_t> rpin;
    nodeID=llabs(nodeID);
    rpin.reserve(reads_in_node[nodeID].size()*2);
    for (auto &rm:reads_in_node[nodeID]){
        rpin.emplace_back(rm.read_id);
        auto other=(rm.read_id%2==1 ? rm.read_id+1:rm.read_id-1);
        if (read_to_node[other]!=nodeID) rpin.emplace_back(other); //only add the other side if it is not in the node.
    }
    std::sort(rpin.begin(),rpin.end());
    return rpin;
}