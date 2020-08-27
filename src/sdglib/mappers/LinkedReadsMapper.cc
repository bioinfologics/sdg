//
// Created by Bernardo Clavijo (EI) on 03/11/2017.
//
#include "LinkedReadsMapper.hpp"
#include <iostream>
#include <iomanip>
#include <cassert>
#include <atomic>
#include <sstream>
#include "sdglib/factories/KMerIDXFactory.hpp"
#include "sdglib/mappers/PerfectMatcher.hpp"
#include <sdglib/datastores/LinkedReadsDatastore.hpp>
#include <sdglib/utilities/omp_safe.hpp>
#include <sdglib/workspace/WorkSpace.hpp>
#include <sdglib/utilities/io_helpers.hpp>
#include <sdglib/datastores/ReadSequenceBuffer.hpp>

const sdgVersion_t LinkedReadsMapper::min_compat = 0x0003;

LinkedReadsMapper::LinkedReadsMapper(WorkSpace &_ws, LinkedReadsDatastore &_datastore) :
ws(_ws),
datastore(_datastore)
{
    reads_in_node.resize(ws.sdg.nodes.size());
}
std::string LinkedReadsMapper::ls(int level,bool recursive) const {
    std::stringstream ss;
    std::string spacer(2 * level, ' ');
    uint64_t mapped=0,unmapped=0;
    for (auto &rtn:read_to_node) if (rtn!=0) ++mapped; else ++unmapped;
    --unmapped;//discard read 0
    ss << spacer << "Linked Reads Mapper: "<<mapped<<" mapped reads, "<<unmapped<<" unmapped" << std::endl;
    return ss.str();
}

void LinkedReadsMapper::write(std::ofstream &output_file) {
    //read-to-node
    uint64_t count=read_to_node.size();
    output_file.write((char *) &SDG_MAGIC, sizeof(SDG_MAGIC));
    output_file.write((char *) &SDG_VN, sizeof(SDG_VN));
    SDG_FILETYPE type(LinkedMap_FT);
    output_file.write((char *) &type, sizeof(type));

    sdglib::write_flat_vector(output_file, read_to_node);
    //mappings
    sdglib::write_flat_vectorvector(output_file, reads_in_node);

    sdglib::write_flat_vectorvector(output_file, tag_neighbours);
}

void LinkedReadsMapper::read(std::ifstream &input_file) {
    uint64_t count;

    sdgMagic_t magic;
    sdgVersion_t version;
    SDG_FILETYPE type;
    input_file.read((char *) &magic, sizeof(magic));
    input_file.read((char *) &version, sizeof(version));
    input_file.read((char *) &type, sizeof(type));

    if (magic != SDG_MAGIC) {
        throw std::runtime_error("Magic number not present in the LinkedReadMap file");
    }

    if (version < min_compat) {
        throw std::runtime_error("LinkedReadMap file version: " + std::to_string(version) + " is not compatible with " + std::to_string(min_compat));
    }

    if (type != LinkedMap_FT) {
        throw std::runtime_error("File type supplied: " + std::to_string(type) + " is not compatible with LinkedMap_FT");
    }


    sdglib::read_flat_vector(input_file, read_to_node);

    sdglib::read_flat_vectorvector(input_file, reads_in_node);

    sdglib::read_flat_vectorvector(input_file, tag_neighbours);

    populate_orientation();


}

void LinkedReadsMapper::dump_readpaths(std::string filename) {
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

void LinkedReadsMapper::load_readpaths(std::string filename) {
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

void LinkedReadsMapper::remap_all_reads() {
    for (auto &rtn:read_to_node) rtn=0;
    for (auto &rin:reads_in_node) rin.clear();
    map_reads();
}

void LinkedReadsMapper::map_reads(const std::unordered_set<uint64_t> &reads_to_remap) {
    UniqueKmerIndex ukindex(ws.sdg);
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
        auto blrs=ReadSequenceBuffer(datastore,128*1024,260);
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

                for (auto &rk:readkmers) {
                    auto nk = ukindex.find(rk.kmer);
                    if (ukindex.end()!=nk) {
                        //get the node just as node
                        sgNodeID_t nknode = llabs(nk->second.node); // nk->second is the graphStrandPosition node is the node id of that
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
    sdglib::OutputLog(sdglib::LogLevels::INFO)<<"Reads mapped: "<<mapped_count<<" / "<<total_count<<" ("<<multimap_count<<" multi-mapped)"<<std::endl;
#pragma omp parallel for
    for (sgNodeID_t n=1;n<reads_in_node.size();++n){
        std::sort(reads_in_node[n].begin(),reads_in_node[n].end());
    }
    populate_orientation();
}

void LinkedReadsMapper::remap_all_reads63() {
    for (auto &rtn:read_to_node) rtn=0;
    for (auto &rin:reads_in_node) rin.clear();
    map_reads63();
}

void LinkedReadsMapper::path_reads(uint8_t k,int _filter) {
    //TODO: when a path jumps nodes implied by overlaps included nodes, add them
    //read pather
    // - Starts from a (multi)63-mer match -> add all starts to PME
    // - extends and check if there's a best path.
    // - uses the index to check if there is further hits, starting from the kmer that includes 1bp beyond the extended match. (cycle)
    sgNodeID_t match_node;
    int64_t match_offset;
    read_paths.clear();
    read_paths.resize(datastore.size()+1);
    if (k<=31) {
        NKmerIndex nki(ws.sdg, k, _filter);//TODO: hardcoded parameters!!
        sdglib::OutputLog() << "Index created!" << std::endl;
        //if (last_read==0) last_read=datastore.size();
#pragma omp parallel shared(nki)
        {
            StreamKmerFactory skf(k);
            std::vector<std::pair<bool, uint64_t >> read_kmers; //TODO: type should be defined in the kmerisers
            std::vector<std::vector<std::pair<int32_t, int32_t >>> kmer_matches; //TODO: type should be defined in the index class
            std::vector<sgNodeID_t> rp;
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
                pme.set_read(seq);
                for (auto rki = 0; rki < read_kmers.size(); ++rki) {
                    auto kmatch = nki.find(read_kmers[rki].second);
                    if (kmatch != nki.end() and kmatch->kmer == read_kmers[rki].second) {
                        pme.reset();
//                    std::cout<<std::endl<<"PME reset done"<<std::endl;
//                    std::cout<<std::endl<<"read kmer is at "<<rki<<" in "<<(read_kmers[rki].first ? "FW":"REV")<<" orientation"<<std::endl;
                        for (; kmatch != nki.end() and kmatch->kmer == read_kmers[rki].second; ++kmatch) {
//                        std::cout<<" match to "<<kmatch->contigID<<":"<<kmatch->offset<<std::endl;
                            auto contig = kmatch->contigID;
                            int64_t pos = kmatch->offset - 1;
                            if (pos < 0) {
                                contig = -contig;
                                pos = -pos - 2;
                            }
                            if (!read_kmers[rki].first) contig = -contig;
                            pme.add_starting_match(contig, rki, pos);
                        }
                        pme.extend_fw();
                        pme.set_best_path();
                        const auto &pmebp = pme.best_path;
                        if (pme.last_readpos == datastore.readsize - 1 or
                            pme.last_readpos - rki >= 63) {//TODO: hardcoded 63bp hit or end
                            if (rp.empty() and not pmebp.empty()) offset = pme.best_path_offset;
                            for (const auto &nid:pmebp)
                                if (rp.empty() or nid != rp.back()) rp.emplace_back(nid);
                            if (!pmebp.empty())
                                rki = pme.last_readpos + 1 - k; //avoid extra index lookups for kmers already used once
                        }
//                    //TODO: matches shold be extended left to avoid unneeded indeterminaiton when an error occurrs in an overlap region and the new hit matches a further part of the genome.
//                    std::cout<<"rki after match extension: "<<rki<<" / "<<read_kmers.size()<<std::endl;
                    }

                }
                read_paths[rid].path = rp;
                read_paths[rid].offset = offset;
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
                            if (rp.empty() and not pmebp.empty()) offset=pme.best_path_offset;
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

void LinkedReadsMapper::map_reads63(const std::unordered_set<uint64_t> &reads_to_remap) {
    Unique63merIndex ukindex(ws.sdg);
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
        auto blrs=ReadSequenceBuffer(datastore,128*1024,260);
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
    sdglib::OutputLog(sdglib::LogLevels::INFO)<<"Reads mapped: "<<mapped_count<<" / "<<total_count<<" ("<<multimap_count<<" multi-mapped)"<<std::endl;
#pragma omp parallel for
    for (sgNodeID_t n=1;n<reads_in_node.size();++n){
        std::sort(reads_in_node[n].begin(),reads_in_node[n].end());
    }
    populate_orientation();
}

void LinkedReadsMapper::print_status() const {
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

std::set<LinkedTag> LinkedReadsMapper::get_node_tags(sgNodeID_t n) {
    std::set<LinkedTag> tags;
    for (auto &rm:reads_in_node[(n>0?n:-n)])
        tags.insert(datastore.get_read_tag(rm.read_id));
    if (tags.count(0)>0) tags.erase(0);
    return tags;
}

std::map<LinkedTag, std::vector<sgNodeID_t>> LinkedReadsMapper::get_tag_nodes(uint32_t min_nodes,
                                                                             const std::vector<bool> &selected_nodes) {
    //Approach: node->tags->nodes (checks how many different tags join this node to every other node).

    std::map<LinkedTag, std::vector<sgNodeID_t>> tag_nodes;
    LinkedTag curr_tag=0,new_tag=0;
    std::set<sgNodeID_t> curr_nodes;
    sgNodeID_t curr_node;
    for (size_t i=1;i<read_to_node.size();++i){
        new_tag=datastore.get_read_tag(i);
        curr_node=llabs(read_to_node[i]);
        if (new_tag==curr_tag){
           if (curr_node!=0 and (selected_nodes.empty() or selected_nodes[curr_node]))
               curr_nodes.insert(curr_node);
        }
        else {
            if (curr_nodes.size()>=min_nodes and curr_tag>0){
                tag_nodes[curr_tag].reserve(curr_nodes.size());
                for (auto &n:curr_nodes)
                    tag_nodes[curr_tag].emplace_back(n);
            }
            curr_tag=new_tag;
            curr_nodes.clear();
            if (curr_node!=0 and (selected_nodes.empty() or selected_nodes[curr_node]))
                curr_nodes.insert(curr_node);
        }
    }
    if (curr_nodes.size()>=min_nodes and curr_tag>0){
        tag_nodes[curr_tag].reserve(curr_nodes.size());
        for (auto &n:curr_nodes)
            tag_nodes[curr_tag].emplace_back(n);
    }

    return tag_nodes;
}

std::vector<std::pair<sgNodeID_t , sgNodeID_t >> LinkedReadsMapper::get_tag_neighbour_nodes(uint32_t min_shared,
                                                                                           const std::vector<bool> &selected_nodes) {
    std::vector<std::pair<sgNodeID_t , sgNodeID_t >> tns;

    auto nodes_in_tags= get_tag_nodes(2, selected_nodes);

    std::unordered_map<sgNodeID_t , std::vector<LinkedTag>> tags_in_node;
    sdglib::OutputLog()<<"There are "<<nodes_in_tags.size()<<" informative tags"<<std::endl;
    for (auto &t:nodes_in_tags){
        for (auto n:t.second) tags_in_node[n].push_back(t.first);
    }

    sdglib::OutputLog()<<"all structures ready"<<std::endl;
#pragma omp parallel firstprivate(tags_in_node,nodes_in_tags)
    {
        std::vector<uint32_t> shared_with(ws.sdg.nodes.size());
        std::vector<std::pair<sgNodeID_t, sgNodeID_t >> tnsl;
#pragma omp for schedule(static,100)
        for (auto n = 1; n < ws.sdg.nodes.size(); ++n) {
            if (!selected_nodes[n]) continue;
            if (tags_in_node[n].size()<min_shared) continue;

            bzero(shared_with.data(), sizeof(uint32_t) * shared_with.size()); //clearing with bzero is the fastest way?
            for (auto t:tags_in_node[n])
                for (auto n:nodes_in_tags[t]) ++shared_with[n];
            for (auto i = n + 1; i < shared_with.size(); ++i) {
                if (shared_with[i] >= min_shared) tnsl.emplace_back(n, i);
            }

        }
#pragma omp critical
        tns.insert(tns.end(),tnsl.begin(),tnsl.end());
    }
    return tns;
}

void LinkedReadsMapper::compute_all_tag_neighbours(int min_size, float min_score, int min_mapped_reads_per_tag) {
    //TODO: parallel for this needs to be done by tag, but tag start-end needs to be computed before

    //first - create start-end for tags that have > X reads
    std::vector<std::pair<uint64_t,uint64_t>> tag_start_end;
    LinkedTag current_tag=0;
    uint64_t mapped_reads=0, tag_firstpost=0;
    for (size_t i=1;i<read_to_node.size();++i){
        if (datastore.get_read_tag(i)!=current_tag){
            if (current_tag!=0 and mapped_reads>=min_mapped_reads_per_tag){
                tag_start_end.emplace_back(tag_firstpost,i);
            }
            mapped_reads=0;
            current_tag=datastore.get_read_tag(i);
            tag_firstpost=i;
        }
        if (read_to_node[i]!=0) ++mapped_reads;
    }
    if (current_tag!=0 and mapped_reads>=min_mapped_reads_per_tag){
        tag_start_end.emplace_back(tag_firstpost,read_to_node.size());
    }

    //scores[source][dest]= count of reads of source which have tags also present in dest
    std::vector<std::unordered_map<sgNodeID_t,uint32_t>> scores(ws.sdg.nodes.size());
    sdglib::OutputLog()<<"Counting into individual maps (parallel)..."<<std::endl;
    {
        //local tally
        std::vector<std::unordered_map<sgNodeID_t,uint32_t>> local_scores(ws.sdg.nodes.size());
        //map with counts of how many times the tag appears on every node
        std::unordered_map<sgNodeID_t, uint32_t > node_readcount; //is it not easier to use a vector? - it is sparse
        for(auto tidx=0;tidx<tag_start_end.size();++tidx){ //for each tag's reads
            for(auto i=tag_start_end[tidx].first;i<tag_start_end[tidx].second;++i) {
                if (read_to_node[i] != 0) {
                    ++node_readcount[llabs(read_to_node[i])];
                }
            }
            if (node_readcount.size()<500) { //tag has reads in less than 500 nodes
                for (auto &n1:node_readcount) {
                    if (ws.sdg.nodes[n1.first].sequence.size() >= min_size) {
                        for (auto &n2:node_readcount) {
                            if (ws.sdg.nodes[n2.first].sequence.size() >= min_size) {
                                scores[n1.first][n2.first] += n1.second; //add the number of reads shared in this tag to the total score
                            }
                        }
                    }
                }
            }
            node_readcount.clear();
        }
        {
            for (auto n=0;n<local_scores.size();++n)
                for (auto ns:local_scores[n])
                    scores[n][ns.first]+=ns.second;
        }
    }

    sdglib::OutputLog()<<"... copying to result vector..."<<std::endl;
    //now flatten the map into its first dimension.
    tag_neighbours.clear();
    tag_neighbours.resize(ws.sdg.nodes.size());
    for (auto i=1;i<ws.sdg.nodes.size();++i){
        for (auto &s:scores[i]) {
            float tag_score = ((float) s.second) / scores[i][i];
            if (tag_score >= min_score)
                tag_neighbours[i].emplace_back(s.first, tag_score);
        }
    }
    sdglib::OutputLog()<<"...DONE!"<<std::endl;
}

LinkedReadsMapper& LinkedReadsMapper::operator=(const LinkedReadsMapper &other) {
    if (&ws.sdg != &other.ws.sdg and &datastore != &other.datastore) { throw std::runtime_error("Can only LinkedReadsMapper from the same SequenceDistanceGraph and LinkedReadsDatastore"); }
    if (&other == this) {
        return *this;
    }
    reads_in_node = other.reads_in_node;
    read_to_node = other.read_to_node;
    return *this;
}

void LinkedReadsMapper::write_tag_neighbours(std::string filename) {
    std::ofstream output_file(filename.c_str());
    uint64_t count;
    count =tag_neighbours.size();
    output_file.write((const char *) &count,sizeof(count));
    for (auto i=0;i<count;++i) {
        uint64_t mcount=tag_neighbours[i].size();
        output_file.write((const char *) &mcount,sizeof(mcount));
        output_file.write((const char *) tag_neighbours[i].data(), sizeof(TagNeighbour) * mcount);
    }
}

void LinkedReadsMapper::read_tag_neighbours(std::string filename) {
    std::ifstream input_file(filename.c_str());
    uint64_t count;
    input_file.read(( char *) &count,sizeof(count));
    tag_neighbours.resize(count);
    for (auto i=0;i<count;++i) {
        uint64_t mcount;
        input_file.read(( char *) &mcount,sizeof(mcount));
        tag_neighbours[i].resize(mcount);
        input_file.read(( char *) tag_neighbours[i].data(), sizeof(TagNeighbour) * mcount);
    }
}

std::ostream &operator<<(std::ostream &os, const LinkedReadsMapper &lirm) {
    uint64_t mapped=0,unmapped=0;
    for (auto &rtn:lirm.read_to_node) if (rtn!=0) ++mapped; else ++unmapped;
    --unmapped;//discard read 0
    os << "Linked Reads Mapper: "<< mapped<<" mapped reads, "<< unmapped<<" unmapped";
    return os;
}

void LinkedReadsMapper::populate_orientation() {
    read_direction_in_node.clear();
    read_direction_in_node.resize(read_to_node.size());
    for (auto & nreads:reads_in_node){
        for (auto &rm:nreads){
            read_direction_in_node[rm.read_id]=rm.rev;
        }
    }
}
