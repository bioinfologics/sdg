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
    output_file.write((char *) &SDG_MAGIC, sizeof(SDG_MAGIC));
    output_file.write((char *) &SDG_VN, sizeof(SDG_VN));
    SDG_FILETYPE type(PairedMap_FT2);
    output_file.write((char *) &type, sizeof(type));
    dump_readpaths(output_file);
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

    if (type == PairedMap_FT) {
        sdglib::OutputLog()<<"WARNING: reading old workspace without PE paths"<<std::endl;
        sdglib::read_flat_vector(input_file, read_to_node);
        sdglib::read_flat_vectorvector(input_file, reads_in_node);
    }
    else if(type == PairedMap_FT2){
        load_readpaths(input_file);
    }
    else throw std::runtime_error("PairedReadsMapper file Incompatible file type");
}


void PairedReadsMapper::dump_readpaths(std::ofstream &opf) {
    uint64_t count;

    //dump operations
    count = read_paths.size();
    opf.write((char *) &count, sizeof(count));
    for (const auto &p:read_paths) {
        opf.write((char *) &p.offset, sizeof(p.offset));
        sdglib::write_flat_vector(opf, p.path);
    }
    sdglib::write_flat_vectorvector(opf,paths_in_node);
}

void PairedReadsMapper::load_readpaths(std::ifstream & ipf) {
    uint64_t count;

    //dump operations
    ipf.read((char *) &count, sizeof(count));
    read_paths.resize(count);
    for (auto &p:read_paths) {
        ipf.read((char *) &p.offset, sizeof(p.offset));
        auto &pp=p.path;
        sdglib::read_flat_vector(ipf,pp);
    }
    sdglib::read_flat_vectorvector(ipf,paths_in_node);
}

void PairedReadsMapper::path_reads(uint8_t k,int _filter,bool fill_offsets) {
    //TODO: when a path jumps nodes implied by overlaps included nodes, add them
    //read pather
    // - Starts from a (multi)63-mer match -> add all starts to PME
    // - extends and check if there's a best path.
    // - uses the index to check if there is further hits, starting from the kmer that includes 1bp beyond the extended match. (cycle)
    sgNodeID_t match_node;
    int64_t match_offset;
    read_paths.clear();
    read_paths.resize(datastore.size()+1);
    if (fill_offsets) read_path_offsets.resize(datastore.size()+1);
    if (k<=31) {
        NKmerIndex nki(ws.sdg, k, _filter);//TODO: hardcoded parameters!!
        sdglib::OutputLog() << "64bit index created!" << std::endl;
        //if (last_read==0) last_read=datastore.size();
#pragma omp parallel shared(nki)
        {
            StreamKmerFactory skf(k);
            std::vector<std::pair<bool, uint64_t >> read_kmers; //TODO: type should be defined in the kmerisers
            std::vector<std::vector<std::pair<int32_t, int32_t >>> kmer_matches; //TODO: type should be defined in the index class
            std::vector<sgNodeID_t> rp;
            std::vector<std::pair<uint32_t,uint32_t>> ro;
            PerfectMatchExtender pme(ws.sdg, k);
            ReadSequenceBuffer sequence_reader(datastore);
            uint32_t offset;
#pragma omp for
            for (uint64_t rid = 1; rid <= datastore.size(); ++rid) {
//        for (auto rid=937;rid<=937;++rid){
                //sdt::cout<<std::endl<<"Processing read "<<rid<<std::endl;
                offset = 0;
                const auto seq = std::string(sequence_reader.get_read_sequence(rid));
                read_kmers.clear();
                skf.produce_all_kmers(seq.c_str(), read_kmers);
                rp.clear();
                ro.clear();
                pme.set_read(seq);
                for (auto rki = 0; rki < read_kmers.size(); ++rki) {
                    auto kmatch = nki.find(read_kmers[rki].second);
                    if (kmatch != nki.end() and kmatch->kmer == read_kmers[rki].second) {
                        pme.reset();
//                    std::cout<<std::endl<<"PME reset done"<<std::endl;
//                    std::cout<<std::endl<<"read kmer is at "<<rki<<" in "<<(read_kmers[rki].first ? "FW":"REV")<<" orientation"<<std::endl;
                        for (; kmatch != nki.end() and kmatch->kmer == read_kmers[rki].second; ++kmatch) {
//                            std::cout<<" match to "<<kmatch->contigID<<":"<<kmatch->offset<<std::endl;
                            auto contig = kmatch->contigID;
                            int64_t pos = kmatch->offset - 1;
                            if (pos < 0) {
                                contig = -contig;
                                pos = -pos - 2;
                            }
                            if (!read_kmers[rki].first) contig = -contig;
//                            std::cout<<"PME starting "<<rki<<" -> "<<contig<<":"<<pos<<std::endl;
                            pme.add_starting_match(contig, rki, pos);
                        }
                        pme.extend_fw();
                        if (fill_offsets) pme.set_best_path(true);
                        else pme.set_best_path();
                        const auto &pmebp = pme.best_path;
//                        std::cout<<"Best path:";
//                        for (auto n:pmebp) std::cout<<" "<<n;
//                        std::cout<<std::endl;
//                        std::cout<<"last read pos: "<<pme.last_readpos<<std::endl;
                        if (not pmebp.empty()) {
                            if (rp.empty() and not pmebp.empty()) offset = pme.best_path_offset;
                            for (uint64_t pnidx=0;pnidx<pmebp.size();++pnidx) {
                                const auto &nid=pmebp[pnidx];
                                if (rp.empty() or nid != rp.back()) {
                                    rp.emplace_back(nid);
                                    if (fill_offsets) ro.emplace_back(pme.best_path_offsets[pnidx]);
                                }
                            }
                            if (!pmebp.empty())
                                rki = pme.last_readpos + 1 - k; //avoid extra index lookups for kmers already used once
                        }
//                    //TODO: matches shold be extended left to avoid unneeded indeterminaiton when an error occurrs in an overlap region and the new hit matches a further part of the genome.
//                    std::cout<<"rki after match extension: "<<rki<<" / "<<read_kmers.size()<<std::endl;
                    }

                }
                read_paths[rid].path = rp;
                read_paths[rid].offset = offset;
                if (fill_offsets) read_path_offsets[rid]=ro;
                //TODO: keep the offset of the first match!
            }
        }
    }
    else if (k<64) {
        NKmerIndex128 nki(ws.sdg,k,_filter);//TODO: hardcoded parameters!!
        sdglib::OutputLog()<<"Index created!"<<std::endl;
        //if (last_read==0) last_read=datastore.size();
#pragma omp parallel shared(nki)
        {
            StreamKmerFactory128 skf(k);
            std::vector<std::pair<bool,__uint128_t >> read_kmers; //TODO: type should be defined in the kmerisers
            std::vector<std::vector<std::pair<int32_t,int32_t >>> kmer_matches; //TODO: type should be defined in the index class
            std::vector<sgNodeID_t> rp;
            std::vector<std::pair<uint32_t,uint32_t>> ro;
            PerfectMatchExtender pme(ws.sdg,k);
            ReadSequenceBuffer sequence_reader(datastore);
            uint32_t offset;
#pragma omp for
            for (uint64_t rid=1;rid<=datastore.size();++rid){
//        for (auto rid=937;rid<=937;++rid){
                //sdt::cout<<std::endl<<"Processing read "<<rid<<std::endl;
                offset=0;
                const auto seq=std::string(sequence_reader.get_read_sequence(rid));
                read_kmers.clear();
                skf.produce_all_kmers(seq.c_str(),read_kmers);
                rp.clear();
                ro.clear();
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
                        if (fill_offsets) pme.set_best_path(true);
                        else pme.set_best_path();
                        const auto &pmebp = pme.best_path;
                        if (not pmebp.empty()) {
                            if (rp.empty() and not pmebp.empty()) offset = pme.best_path_offset;
                            for (uint64_t pnidx=0;pnidx<pmebp.size();++pnidx) {
                                const auto &nid=pmebp[pnidx];
                                if (rp.empty() or nid != rp.back()) {
                                    rp.emplace_back(nid);
                                    if (fill_offsets) ro.emplace_back(pme.best_path_offsets[pnidx]);
                                }
                            }
                            if (!pmebp.empty())
                                rki = pme.last_readpos + 1 - k; //avoid extra index lookups for kmers already used once
                        }
//                    //TODO: matches shold be extended left to avoid unneeded indeterminaiton when an error occurrs in an overlap region and the new hit matches a further part of the genome.
//                    std::cout<<"rki after match extension: "<<rki<<" / "<<read_kmers.size()<<std::endl;
                    }

                }
                read_paths[rid].path=rp;
                read_paths[rid].offset=offset;
                //TODO: keep the offset of the first match!
            }
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

}

bool inline are_complement(const char &A,const char &B){
    return (A=='A' and B=='T') or (A=='C' and B=='G') or (A=='G' and B=='C') or (A=='T' and B=='A');
}

std::vector<sgNodeID_t> PairedReadsMapper::path_fw(seqID_t read_id, sgNodeID_t node, bool use_pair, bool collapse_pair) const{
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

std::vector<std::vector<sgNodeID_t> > PairedReadsMapper::all_paths_fw(sgNodeID_t node, bool use_pair, bool collapse_pair) const {
    if (llabs(node)>paths_in_node.size()) return {};
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

std::ostream &operator<<(std::ostream &os, const PairedReadsMapper &prm) {
    uint64_t mapped=0,unmapped=0;
    for (auto &rtn:prm.read_to_node) if (rtn!=0) ++mapped; else ++unmapped;
    if(unmapped>0)--unmapped;//discard read 0

    os << "Paired Reads Mapper: "<<mapped<<" mapped reads, "<<unmapped<<" unmapped";
    return os;
}

sgNodeID_t PairedReadsMapper::get_node_inmediate_neighbours(sgNodeID_t node){
    // 0 returning node is equivalent to None (in python)
    sgNodeID_t next_node=0;
    std::unordered_map<sgNodeID_t, int> counter;
    for(const auto &path: all_paths_fw(node)){
        for (const auto &node: path){
            counter[node]++;
        }
    }

    std::vector<std::pair<sgNodeID_t, int>> count_vector;
    for (const auto &x: counter){
        count_vector.push_back(x);
    }
    std::sort(count_vector.begin(), count_vector.end(), [](std::pair<sgNodeID_t, int>& a, std::pair<sgNodeID_t, int>& b){return a.second>b.second;});

    for (const auto& cand: count_vector){
        if (cand.first != 0 and cand.second>10){
            return cand.first;
        }
    }
    return next_node;
}

std::vector<int64_t> PairedReadsMapper::get_paths_in_node(sgNodeID_t nid){
    if (nid>=0) return paths_in_node[nid];
    std::vector<int64_t> p;
    p.reserve(paths_in_node[-nid].size());
    for (auto &pin:paths_in_node[-nid]) p.emplace_back(-pin);
    return p;
}