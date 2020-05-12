//
// Created by Bernardo Clavijo (EI) on 03/11/2017.
//
#include "PairedReadsMapper.hpp"
#include <iomanip>
#include <cassert>
#include <atomic>
#include <sdglib/utilities/omp_safe.hpp>
#include <fstream>
#include <iostream>
#include <sstream>
#include <sdglib/workspace/WorkSpace.hpp>
#include <sdglib/utilities/io_helpers.hpp>
#include <sdglib/views/NodeView.hpp>
#include <sdglib/mappers/PerfectMatcher.hpp>

const sdgVersion_t PairedReadsMapper::min_compat = 0x0003;

PairedReadsMapper::PairedReadsMapper(WorkSpace &_ws, PairedReadsDatastore &_datastore) :
        ws(_ws),
        datastore(_datastore)
{
    reads_in_node.resize(ws.sdg.nodes.size());
}

std::string PairedReadsMapper::ls(int level,bool recursive) const {
    std::stringstream ss;
    std::string spacer(2 * level, ' ');
    uint64_t mapped=0,unmapped=0;
    for (auto &rtn:read_to_node) if (rtn!=0) ++mapped; else ++unmapped;
    if(unmapped>0)--unmapped;//discard read 0
    ss << spacer << "Paired Reads Mapper: "<<mapped<<" mapped reads, "<<unmapped<<" unmapped" << std::endl;
    return ss.str();
}

void PairedReadsMapper::write(std::ofstream &output_file) {
    //read-to-node
    output_file.write((char *) &SDG_MAGIC, sizeof(SDG_MAGIC));
    output_file.write((char *) &SDG_VN, sizeof(SDG_VN));
    SDG_FILETYPE type(PairedMap_FT);
    output_file.write((char *) &type, sizeof(type));

    sdglib::write_flat_vector(output_file, read_to_node);
    //mappings
    sdglib::write_flat_vectorvector(output_file, reads_in_node);
}

void PairedReadsMapper::read(std::ifstream &input_file) {
    sdgMagic_t magic;
    sdgVersion_t version;
    SDG_FILETYPE type;
    input_file.read((char *) &magic, sizeof(magic));
    input_file.read((char *) &version, sizeof(version));
    input_file.read((char *) &type, sizeof(type));

    if (magic != SDG_MAGIC) {
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


void PairedReadsMapper::dump_readpaths(std::string filename) {
    std::ofstream opf(filename);
    uint64_t count;

    //dump operations
    count = read_paths.size();
    opf.write((char *) &count, sizeof(count));
    for (const auto &p:read_paths) {
        opf.write((char *) &p.offset, sizeof(p.offset));
        sdglib::write_flat_vector(opf, p.path);
    }
}

void PairedReadsMapper::load_readpaths(std::string filename) {
    std::ifstream ipf(filename);
    uint64_t count;

    //dump operations
    ipf.read((char *) &count, sizeof(count));
    read_paths.resize(count);
    for (auto &p:read_paths) {
        ipf.read((char *) &p.offset, sizeof(p.offset));
        auto &pp=p.path;
        sdglib::read_flat_vector(ipf,pp);
    }

    sdglib::OutputLog(sdglib::LogLevels::INFO)<<"Reserving paths_in_node"<<std::endl;
    //first pass, size computing for vector allocation;
    std::vector<uint64_t> path_counts(ws.sdg.nodes.size());
    for (auto i=0;i<read_paths.size();++i){
        const auto &p=read_paths[i];
        for (const auto &n:p.path){
            ++path_counts[llabs(n)];
        }
    }

    paths_in_node.clear();
    paths_in_node.resize(path_counts.size());
    for (auto i=0;i<path_counts.size();++i) paths_in_node[i].reserve(path_counts[i]);
    sdglib::OutputLog(sdglib::LogLevels::INFO)<<"Filling paths_in_node"<<std::endl;
    for (auto i=0;i<read_paths.size();++i){
        const auto &p=read_paths[i];
        for (const auto &n:p.path){
            auto pid=(n<0 ? -i : i);
            if (paths_in_node[llabs(n)].empty() or paths_in_node[llabs(n)].back()!=pid) paths_in_node[llabs(n)].emplace_back(pid);
        }
    }
    sdglib::OutputLog(sdglib::LogLevels::INFO)<<"paths_in_node filled"<<std::endl;
}

void PairedReadsMapper::remap_all_reads() {
    for (auto &rtn:read_to_node) rtn=0;
    for (auto &rin:reads_in_node) rin.clear();
    map_reads();
}

void PairedReadsMapper::map_reads(const std::unordered_set<uint64_t> &reads_to_remap) {
    UniqueKmerIndex ukindex(ws.sdg);
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
        auto blrs=ReadSequenceBuffer(datastore,128*1024,datastore.readsize*2+2);
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
                    auto nk = ukindex.find(rk.kmer);
                    if (ukindex.end()!=nk) {
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
    Unique63merIndex ukindex(ws.sdg);
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
        auto blrs=ReadSequenceBuffer(datastore,128*1024,datastore.readsize*2+2);
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
                    auto nk = ukindex.find(rk.kmer);
                    if (ukindex.end()!=nk) {
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



void PairedReadsMapper::path_reads() {
    //TODO: when a path jumps nodes implied by overlaps included nodes, add them
    //read pather
    // - Starts from a (multi)63-mer match -> add all starts to PME
    // - extends and check if there's a best path.
    // - uses the index to check if there is further hits, starting from the kmer that includes 1bp beyond the extended match. (cycle)
    sgNodeID_t match_node;
    int64_t match_offset;
    read_paths.clear();
    read_paths.resize(datastore.size()+1);
    uint8_t k=31;
    NKmerIndex nki(ws.sdg,k,50);//TODO: hardcoded parameters!!
    //sdt::cout<<"Index created!"<<std::endl;
    //if (last_read==0) last_read=datastore.size();
#pragma omp parallel shared(nki)
    {
        StreamKmerFactory skf(k);
        std::vector<std::pair<bool,uint64_t >> read_kmers; //TODO: type should be defined in the kmerisers
        std::vector<std::vector<std::pair<int32_t,int32_t >>> kmer_matches; //TODO: type should be defined in the index class
        std::vector<sgNodeID_t> rp;
        PerfectMatchExtender pme(ws.sdg,k);
        ReadSequenceBuffer sequence_reader(datastore);
#pragma omp for
        for (auto rid=1;rid<=datastore.size();++rid){
//        for (auto rid=937;rid<=937;++rid){
            //sdt::cout<<std::endl<<"Processing read "<<rid<<std::endl;
            const auto seq=std::string(sequence_reader.get_read_sequence(rid));
            read_kmers.clear();
            skf.produce_all_kmers(seq.c_str(),read_kmers);
            rp.clear();
            pme.set_read(seq);
            for (auto rki=0;rki<read_kmers.size();++rki){
                auto kmatch=nki.find(read_kmers[rki].second);
                if (kmatch!=nki.end() and kmatch->kmer==read_kmers[rki].second){
                    pme.reset();
//                    std::cout<<std::endl<<"PME reset done"<<std::endl;
//                    std::cout<<std::endl<<"read kmer is at "<<rki<<" in "<<(read_kmers[rki].first ? "FW":"REV")<<" orientation"<<std::endl;
                    for (;kmatch!=nki.end() and kmatch->kmer==read_kmers[rki].second;++kmatch) {
//                        std::cout<<" match to "<<kmatch->contigID<<":"<<kmatch->offset<<std::endl;
                        auto contig=kmatch->contigID;
                        int64_t pos=kmatch->offset-1;
                        if (pos<0) {
                            contig=-contig;
                            pos=-pos-2;
                        }
                        if (!read_kmers[rki].first) contig=-contig;
                        pme.add_starting_match(contig, rki, pos);
                    }
                    pme.extend_fw();
                    pme.set_best_path();
                    const auto & pmebp=pme.best_path;
                    if (pme.last_readpos==datastore.readsize-1 or pme.last_readpos-rki>=63) {//TODO: hardcoded 63bp hit or end
                        for (const auto &nid:pmebp)
                            if (rp.empty() or nid != rp.back()) rp.emplace_back(nid);
                        if (!pmebp.empty())
                            rki = pme.last_readpos + 1 - k; //avoid extra index lookups for kmers already used once
                    }
//                    //TODO: matches shold be extended left to avoid unneeded indeterminaiton when an error occurrs in an overlap region and the new hit matches a further part of the genome.
//                    std::cout<<"rki after match extension: "<<rki<<" / "<<read_kmers.size()<<std::endl;
                }

            }
            read_paths[rid].path=rp;
            //TODO: keep the offset of the first match!
        }
    }
    sdglib::OutputLog(sdglib::LogLevels::INFO)<<"Creating paths_in_node"<<std::endl;
    //first pass, size computing for vector allocation;
    std::vector<uint64_t> path_counts(ws.sdg.nodes.size());
    for (auto i=0;i<read_paths.size();++i){
        const auto &p=read_paths[i];
        for (const auto &n:p.path){
            ++path_counts[llabs(n)];
        }
    }

    paths_in_node.clear();
    paths_in_node.resize(path_counts.size());
    for (auto i=0;i<path_counts.size();++i) paths_in_node[i].reserve(path_counts[i]);
    for (auto i=0;i<read_paths.size();++i){
        const auto &p=read_paths[i];
        for (const auto &n:p.path){
            auto pid=(n<0 ? -i : i);
            if (paths_in_node[llabs(n)].empty() or paths_in_node[llabs(n)].back()!=pid) paths_in_node[llabs(n)].emplace_back(pid);
        }
    }
    sdglib::OutputLog(sdglib::LogLevels::INFO)<<"pathing finished"<<std::endl;

};

bool inline are_complement(const char &A,const char &B){
    return (A=='A' and B=='T') or (A=='C' and B=='G') or (A=='G' and B=='C') or (A=='T' and B=='A');
}

void PairedReadsMapper::path_reads63() {
    const int k = 63;
    Unique63merIndex ukindex(ws.sdg);
    std::atomic<int64_t> nokmers(0);
    read_paths.clear();
    read_paths.resize(datastore.size()+1);

    /*
     * Read mapping in parallel,
     */
    uint64_t thread_mapped_count[omp_get_max_threads()],thread_total_count[omp_get_max_threads()],thread_multimap_count[omp_get_max_threads()];
    std::vector<ReadMapping> thread_mapping_results[omp_get_max_threads()];
    sdglib::OutputLog()<<"Mapping with "<<omp_get_max_threads()<<" threads"<<std::endl;
#pragma omp parallel
    {
        std::vector<KmerIDX128> readkmers;
        StreamKmerIDXFactory128 skf(63);
        auto blrs=ReadSequenceBuffer(datastore,128*1024,datastore.readsize*2+2);
        auto & private_results=thread_mapping_results[omp_get_thread_num()];
        auto & mapped_count=thread_mapped_count[omp_get_thread_num()];
        auto & total_count=thread_total_count[omp_get_thread_num()];
        auto & multimap_count=thread_multimap_count[omp_get_thread_num()];
        mapped_count=0;
        total_count=0;
        multimap_count=0;
#pragma omp for
        for (uint64_t readID=1;readID<read_paths.size();++readID) {
            auto &path=read_paths[readID].path;
            path.clear();
            //get all kmers from read
            auto seq=blrs.get_read_sequence(readID);
            readkmers.clear();
            skf.produce_all_kmers(seq,readkmers);
            if (readkmers.size()==0) {
                ++nokmers;
            }
            //lookup for kmer, but extend on current sequence until missmatch or end
            for (auto rki=0;rki<readkmers.size();++rki) {
                auto &rk=readkmers[rki];
                auto nk = ukindex.find(rk.kmer);
                if (ukindex.end()!=nk) {
                    //get the node just as node
                    sgNodeID_t nknode = nk->second.node;
                    if (rk.contigID<0) nknode=-nknode;
                    if ( nknode>0 ) {
                        auto cpos = nk->second.pos + k;
                        while (rki < readkmers.size()
                               and cpos < ws.sdg.nodes[nknode].sequence.size()
                               and ws.sdg.nodes[nknode].sequence[cpos] == seq[rki + k]) {
                            ++cpos, ++rki;
                        }
                    }
                    else {
                        int32_t cpos = nk->second.pos-1;
                        while (rki < readkmers.size()
                               and cpos >= 0
                               and are_complement(ws.sdg.nodes[-nknode].sequence[cpos],seq[rki + k])) {
                            --cpos, ++rki;
                        }
                    }
                    if (path.empty() or path.back()!=nknode) path.emplace_back(nknode);
                }
            }
            if (not path.empty()) {
                ++mapped_count;
                if (path.size()>1) ++multimap_count;
            }
            auto tc = ++total_count;
            if (tc % 10000000 == 0) sdglib::OutputLog()<< mapped_count << " / " << tc <<" ("<<++multimap_count<<" multi-node paths)"<< std::endl;
        }
    }
    uint64_t mapped_count=0,total_count=0,multimap_count=0;
    for (auto i=0;i<omp_get_max_threads();++i){
        mapped_count+=thread_mapped_count[i];
        total_count+=thread_total_count[i];
        multimap_count+=thread_multimap_count[i];
    }
    sdglib::OutputLog(sdglib::LogLevels::INFO)<<"Reads without k-mers: "<<nokmers<<std::endl;
    sdglib::OutputLog(sdglib::LogLevels::INFO)<<"Reads mapped: "<<mapped_count<<" / "<<total_count<<" ("<<multimap_count<<" multi-node paths)"<<std::endl;

    sdglib::OutputLog(sdglib::LogLevels::INFO)<<"Creating paths_in_node"<<std::endl;
    //first pass, size computing for vector allocation;
    std::vector<uint64_t> path_counts(ws.sdg.nodes.size());
    for (auto i=0;i<read_paths.size();++i){
        const auto &p=read_paths[i];
        for (const auto &n:p.path){
            ++path_counts[llabs(n)];
        }
    }

    paths_in_node.clear();
    paths_in_node.resize(path_counts.size());
    for (auto i=0;i<path_counts.size();++i) paths_in_node[i].reserve(path_counts[i]);
    for (auto i=0;i<read_paths.size();++i){
        const auto &p=read_paths[i];
        for (const auto &n:p.path){
            auto pid=(n<0 ? -i : i);
            if (paths_in_node[llabs(n)].empty() or paths_in_node[llabs(n)].back()!=pid) paths_in_node[llabs(n)].emplace_back(pid);
        }
    }
    sdglib::OutputLog(sdglib::LogLevels::INFO)<<"pathing finished"<<std::endl;
}

std::vector<sgNodeID_t> PairedReadsMapper::path_fw(seqID_t read_id, sgNodeID_t node, bool use_pair, bool collapse_pair) {
    auto pread_id=datastore.get_read_pair(read_id);

    //check for obvious no-path conditions:
    if (read_paths[llabs(read_id)].path.empty()) return {};
    if (read_paths[llabs(read_id)].path.size()==1) {
        if (not use_pair) return {};
        if (read_paths[llabs(pread_id)].path.empty()) return {};
        if (read_paths[llabs(pread_id)].path.size()==1 and read_paths[llabs(pread_id)].path[0]==(pread_id>0 ? -node:node) ) return {};
    }

    //get read paths in the right orientation - leaves r2p empty if not using pairs
    std::vector<sgNodeID_t> r1p,r2p;
    if (read_id > 0) {
        r1p=read_paths[read_id].path;
        if (use_pair) sdglib::reverse_path(r2p,read_paths[-datastore.get_read_pair(read_id)].path);
    }
    else {
        sdglib::reverse_path(r1p,read_paths[-read_id].path); //no read pair, as the pair comes "BEFORE"
    }

    //if collapsing pairs, find max overlap between end of r1p and begining of r2p
    int ovl=0;
    if (collapse_pair){
        for (ovl=std::min(r1p.size(),r2p.size());ovl>0;--ovl){
            auto p1it=r1p.cend()-ovl, p2it=r2p.cbegin();
            while(p1it<r1p.cend() and *p1it==*p2it) ++p1it,++p2it;
            if (p1it==r1p.cend()) break;
        }
    }

    //discards up to node in r1p
    auto r1it=r1p.cbegin();
    while (r1it<r1p.cend() and *r1it!=node) ++r1it;
    if (r1it==r1p.cend()) return {};
    std::vector<sgNodeID_t> path_fw;
    path_fw.insert(path_fw.end(),++r1it,r1p.cend());

    //adds {0} and r2p if there's anything there
    if (ovl<r2p.size()){
        path_fw.emplace_back(0);
        path_fw.insert(path_fw.end(),r2p.cbegin()+ovl,r2p.cend());
    }
    return path_fw;
}

std::vector<std::vector<sgNodeID_t> > PairedReadsMapper::all_paths_fw(sgNodeID_t node, bool use_pair, bool collapse_pair) {
    std::vector<std::vector<sgNodeID_t> > r;
    std::unordered_set<seqID_t> rids;
    for (auto rid:paths_in_node[llabs(node)] ){
        if (node<0) rid=-rid;
        if (rid<0 or rids.count(datastore.get_read_pair(rid))==0) {
            rids.insert(rid);
            r.emplace_back(path_fw(rid,node,use_pair,collapse_pair));
            if (r.back().empty()) r.pop_back();
        }
    }
    return r;
}

void PairedReadsMapper::print_status() const {
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

std::ostream &operator<<(std::ostream &os, const PairedReadsMapper &prm) {
    uint64_t mapped=0,unmapped=0;
    for (auto &rtn:prm.read_to_node) if (rtn!=0) ++mapped; else ++unmapped;
    if(unmapped>0)--unmapped;//discard read 0

    os << "Paired Reads Mapper: "<<mapped<<" mapped reads, "<<unmapped<<" unmapped";
    return os;
}
