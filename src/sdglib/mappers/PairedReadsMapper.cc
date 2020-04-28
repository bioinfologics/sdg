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

bool inline are_complement(const char &A,const char &B){
    return (A=='A' and B=='T') or (A=='C' and B=='G') or (A=='G' and B=='C') or (A=='T' and B=='A');
}

class PerfectMatchPart{
public:
    void extend(const std::string & readseq,const std::string & nodeseq) {//this jsut grows until it can't and sets the flags.
        //while(i<size(read) and j<size(node) and read[i]==node[j]) ++i (consider RC, maybe just write the conditions appropriately?)
        if (node>0){
            while(read_position<readseq.size()-1 and node_position<nodeseq.size()-1 and nodeseq[node_position+1]==readseq[read_position+1]){
                ++node_position;
                ++read_position;
            }
            completed_node = (node_position==nodeseq.size()-1);
            completed_read = (read_position==readseq.size()-1);
        }
        else {
            while(read_position<readseq.size()-1 and node_position>0 and are_complement(nodeseq[node_position-1],readseq[read_position+1])){
                --node_position;
                ++read_position;
            }
            completed_node = (node_position==0);
            completed_read = (read_position==readseq.size()-1);
        }
        extended=true;
    }

    sgNodeID_t node;

    uint64_t offset;//will only be set if first part
    int64_t previous_part;//will only be set if not first part
    uint64_t read_position;//position of the last matched base in read
    uint64_t node_position;//postiion of the last matched base in node (canonical orientation)

    bool completed_node=false;
    bool completed_read=false;
    bool extended=false;
};

class PerfectMatchExtender{
public:
    PerfectMatchExtender(DistanceGraph & _dg):dg(_dg){};

    void reset(std::string _readseq){
        matchparts.clear();
        readseq=_readseq;
        last_readpos=0;
    }

    void add_starting_match(sgNodeID_t node_id, uint64_t read_offset, uint64_t node_offset, uint8_t _k){ //add a new start
        k=_k;
        matchparts.emplace_back();
        auto &mp=matchparts.back();
        mp.previous_part=-1;
        mp.node=node_id;
        mp.extended=false;
        mp.read_position=read_offset+k-1;
        mp.extended=mp.completed_read=(mp.read_position==readseq.size()-1);
        if (node_id>0) {
            mp.node_position=node_offset+k-1;
            mp.completed_node = (mp.node_position==dg.sdg.get_node_size(node_id)-1);
        }
        else {
            mp.node_position=node_offset;
            mp.completed_node = (node_offset==0);
        }
    }

    void extend_fw(){
        //TODO: this could be extended backwards in the original hit to add diferentiation on overlap-transitioned hits after an error.
        std::cout<<"extend_fw called with "<<matchparts.size()<<" starting matchparts"<<std::endl;
        for(uint64_t next=0;next<matchparts.size();++next){
            //extend, if end of node add all nexts as unextended parts.
            std::cout<<"extending matchpart "<<next<<" to node "<<matchparts[next].node<<" with current readpos="<<matchparts[next].read_position<<" and nodepos="<<matchparts[next].node_position<<std::endl;
            matchparts[next].extend(readseq,dg.sdg.nodes[llabs(matchparts[next].node)].sequence);
            std::cout<<" -> readpos="<<matchparts[next].read_position<<(matchparts[next].completed_read ? " (completed)":"")<<", nodepos="<<matchparts[next].node_position<<(matchparts[next].completed_node ? " (completed)":"")<<std::endl;
            if (matchparts[next].completed_node and not matchparts[next].completed_read){
                for (auto l: dg.get_nodeview(matchparts[next].node).next()){
                    if (l.distance()>-k+1) continue;
                    std::cout<<"extending to next node "<<l.node().node_id()<<std::endl;
                    matchparts.emplace_back();//add next node part, pointing to this one as previous.
                    auto &mp=matchparts.back();
                    mp.node=l.node().node_id();
                    mp.read_position=matchparts[next].read_position+1;
                    if (mp.node>0) {
                        mp.node_position=-l.distance();
                    }
                    else {
                        mp.node_position=dg.sdg.get_node_size(mp.node)+l.distance()-1;
                    }
                    mp.extended=false;
                    mp.previous_part=next;
                }
            }
        }
    }
    std::vector<sgNodeID_t> best_path(){ //set extension_size when returning a path that was extended.
        std::vector<sgNodeID_t> bp;
        //find the matches that goes farthest in the read?
        //more than one? -> unspecific section, return nothing. //TODO: keep them in case many of them define a longest path still?
        uint64_t best_length=0;
        int64_t next_index=0;
        for (auto i=0;i<matchparts.size();++i){
            auto &mp=matchparts[i];
            if (mp.read_position==best_length) next_index=-1;
            else if (mp.read_position>best_length) {
                best_length=mp.read_position;
                next_index=i;
            }
        }
        last_nodepos=matchparts[next_index].node_position;

        while (next_index!=-1) {
            bp.emplace_back(matchparts[next_index].node);
            next_index=matchparts[next_index].previous_part;
        }
        std::reverse(bp.begin(),bp.end());
        if (!bp.empty()) last_readpos=best_length;
        else last_nodepos=0;
        return bp;
    }

    DistanceGraph & dg;
    uint8_t k;
    std::vector<PerfectMatchPart> matchparts;
    uint64_t last_readpos;
    uint64_t last_nodepos;
    std::string readseq;

};

void PairedReadsMapper::path_reads() {
    //TODO: compute overlap regions for nodes to not start a match into an overlap nor accept a match fully on an overlap.
    //read pather
    // - Starts from a (multi)63-mer match -> add all starts to PME
    // - extends and check if there's a best path.
    // - uses the index to check if there is further hits, starting from the kmer that includes 1bp beyond the extended match. (cycle)
    sgNodeID_t match_node;
    int64_t match_offset;
    read_paths.clear();
    read_paths.resize(datastore.size());
    uint8_t k=31;
    NKmerIndex nki(ws.sdg,k,50);//TODO: hardcoded parameters!!
    std::cout<<"Index created!"<<std::endl;
    //if (last_read==0) last_read=datastore.size();
//#pragma omp parallel shared(nki)
    {
        StreamKmerFactory skf(k);
        std::vector<std::pair<bool,uint64_t >> read_kmers; //TODO: type should be defined in the kmerisers
        std::vector<std::vector<std::pair<int32_t,int32_t >>> kmer_matches; //TODO: type should be defined in the index class
        std::vector<sgNodeID_t> rp;
        PerfectMatchExtender pme(ws.sdg);
        ReadSequenceBuffer sequence_reader(datastore);
//#pragma omp for
        //for (auto rid=1;rid<=datastore.size();++rid){
        for (auto rid=2;rid<=2;++rid){
            std::cout<<std::endl<<"Processing read "<<rid<<std::endl;
            auto seq=sequence_reader.get_read_sequence(rid);
            read_kmers.clear();
            skf.produce_all_kmers(seq,read_kmers);
            rp.clear();
            for (auto rki=0;rki<read_kmers.size();++rki){
                auto kmatch=nki.find(read_kmers[rki].second);
                if (kmatch!=nki.end() and kmatch->kmer==read_kmers[rki].second){
                    pme.reset(seq);
                    //std::cout<<std::endl<<"PME reset done"<<std::endl;
                    std::cout<<std::endl<<"read kmer is at "<<rki<<" in "<<(read_kmers[rki].first ? "FW":"REV")<<" orientation"<<std::endl;
                    for (;kmatch!=nki.end() and kmatch->kmer==read_kmers[rki].second;++kmatch) {
                        std::cout<<" match to "<<kmatch->contigID<<":"<<kmatch->offset<<std::endl;
                        auto contig=kmatch->contigID;
                        int64_t pos=kmatch->offset-1;
                        if (pos<0) {
                            contig=-contig;
                            pos=-pos-2;
                        }
                        if (!read_kmers[rki].first) contig=-contig;
                        //std::cout<<"pos = "<<ws.sdg.get_node_size(contig)<<" / "<<pos<<std::endl;
                        //avoid starting a node in overlap zone, TODO:hardcoded for overlaps of 63!!!
                        if ( (contig>0 and ws.sdg.get_node_size(contig)-pos<63) or (contig<0 and pos<63-k) )
                            pme.add_starting_match(contig, rki, pos, k);

                    }
                    pme.extend_fw();
                    auto pmebp=pme.best_path();

                    if (pmebp.size()>1 or (pmebp.size()==1 and ( (kmatch->offset>0 and pme.last_nodepos>63) or (kmatch->offset<0 and  ws.sdg.get_node_size(kmatch->contigID)-pme.last_nodepos>63)))) {
                        for (auto &nid:pmebp)
                            if (rp.empty() or nid != rp.back()) rp.emplace_back(nid);
                    }
                    if (!pmebp.empty()) rki=pme.last_readpos+1-k; //this avoids doing extra index looks for kmers that are already used-up.
                    //TODO: matches shold be extended left to avoid unneeded indeterminaiton when an error occurrs in an overlap region and the new hit matches a further part of the genome.
                    std::cout<<"rki after match extension: "<<rki<<" / "<<read_kmers.size()<<std::endl;
                }

            }
            read_paths[rid].path=rp;
            //TODO: keep the offset of the first match!
        }
    }

};

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

std::vector<sgNodeID_t> PairedReadsMapper::path_fw(int64_t read_id, sgNodeID_t node, bool use_pair) {
    std::vector<sgNodeID_t> path_fw;

    if (read_id>0) {
        auto it=read_paths[read_id].path.cbegin();
        while (it!=read_paths[read_id].path.cend() and *it!=node) ++it;
        if (it==read_paths[read_id].path.cend()) return {};
        else ++it;
        while (it!=read_paths[read_id].path.cend()) path_fw.emplace_back(*it++);
        if (use_pair) {
            path_fw.emplace_back(0);
            auto read2_id = (read_id % 2 ? read_id + 1 : read_id - 1);
            auto itr = read_paths[read2_id].path.crbegin();
            for (; itr != read_paths[read2_id].path.crend() and *itr != -node; ++itr);
            if (itr == read_paths[read2_id].path.crend()) itr = read_paths[read2_id].path.crbegin();
            else ++itr;
            while (itr != read_paths[read2_id].path.crend()) path_fw.emplace_back(-*itr++);
        }
    }
    else {
        auto itr=read_paths[-read_id].path.crbegin();
        while (itr!=read_paths[-read_id].path.crend() and *itr!=-node) ++itr;
        if (itr==read_paths[-read_id].path.crend()) return {};
        else ++itr;
        while (itr!=read_paths[-read_id].path.crend()) path_fw.emplace_back(-*itr++);
        //pair is not FW, so it is not used!
    }
    if (not path_fw.empty() and path_fw.back()==0) path_fw.resize(path_fw.size()-1);
    return path_fw;
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

std::ostream &operator<<(std::ostream &os, const PairedReadsMapper &prm) {
    uint64_t mapped=0,unmapped=0;
    for (auto &rtn:prm.read_to_node) if (rtn!=0) ++mapped; else ++unmapped;
    if(unmapped>0)--unmapped;//discard read 0

    os << "Paired Reads Mapper: "<<mapped<<" mapped reads, "<<unmapped<<" unmapped";
    return os;
}
