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
    output_file.write((char *) &SDG_MAGIC, sizeof(SDG_MAGIC));
    output_file.write((char *) &SDG_VN, sizeof(SDG_VN));
    SDG_FILETYPE type(LinkedMap_FT);
    output_file.write((char *) &type, sizeof(type));
    dump_readpaths(output_file);
    sdglib::write_flat_vectorvector(output_file, tag_neighbours);
}

void LinkedReadsMapper::read(std::ifstream &input_file) {

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

    if (type == LinkedMap_FT){
        load_readpaths(input_file);
        sdglib::read_flat_vectorvector(input_file, tag_neighbours);
    }
    else throw std::runtime_error("PairedReadsMapper file Incompatible file type");

}

void LinkedReadsMapper::dump_readpaths(std::ofstream &opf) {
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

void LinkedReadsMapper::load_readpaths(std::ifstream & ipf) {
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

void LinkedReadsMapper::path_reads(uint8_t k,int _filter,bool fill_offsets) {
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
    for (auto &p:paths_in_node[(n>0?n:-n)])
        tags.insert(datastore.get_read_tag(llabs(p)));
    if (tags.count(0)>0) tags.erase(0);
    return tags;
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
