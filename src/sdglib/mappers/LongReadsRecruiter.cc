//
// Created by Bernardo Clavijo (EI) on 16/01/2020.
//

#include "LongReadsRecruiter.hpp"
#include "PerfectMatcher.hpp"
#include <sdglib/indexers/NKmerIndex.hpp>
#include <atomic>

LongReadsRecruiter::LongReadsRecruiter(SequenceDistanceGraph &_sdg, const LongReadsDatastore &datastore,uint8_t _k=25, uint16_t _f=50):
    sdg(_sdg),datastore(datastore),k(_k),f(_f) {
    reset_recruitment();
};


void LongReadsRecruiter::dump(std::string filename) {
    std::ofstream ofs(filename);
    sdglib::write_flat_vectorvector(ofs,node_reads);
    sdglib::write_flat_vectorvector(ofs,read_perfect_matches);
}

void LongReadsRecruiter::load(std::string filename) {
    std::ifstream ifs(filename);
    sdglib::read_flat_vectorvector(ifs,node_reads);
    sdglib::read_flat_vectorvector(ifs,read_perfect_matches);
}

void LongReadsRecruiter::reset_recruitment() {
    node_reads.clear();
    node_reads.resize(sdg.nodes.size());
    read_perfect_matches.clear();
    read_perfect_matches.resize(datastore.size()+1);
}
//TODO: optimise the in-vector class to do a binary search (right now the kmatches are not well-ordered because of RCs)
template<class T>
inline bool in_vector(const std::vector<T> &V,T VAL) { return std::find(V.begin(),V.end(),VAL)!=V.end();}

void LongReadsRecruiter::perfect_mappings(uint16_t seed_size, uint64_t first_read, uint64_t last_read) {
    if (seed_size<k) seed_size=k;
    NKmerIndex nki(sdg,k,f);
    if (last_read==0) last_read=datastore.size();
#pragma omp parallel shared(nki)
    {
        StreamKmerFactory skf(k);
        std::vector<std::pair<bool,uint64_t >> read_kmers; //TODO: type should be defined in the kmerisers
        std::vector<std::vector<std::pair<int32_t,int32_t >>> kmer_matches; //TODO: type should be defined in the index class
        PerfectMatch pmatch(0,0,0);
        ReadSequenceBuffer sequence_reader(datastore);
#pragma omp for
        for (auto rid=first_read;rid<=last_read;++rid){
            std::vector<PerfectMatch> private_read_perfect_matches;
            auto seq=sequence_reader.get_read_sequence(rid);
            read_kmers.clear();
            skf.produce_all_kmers(seq,read_kmers);
            nki.get_all_kmer_matches(kmer_matches,read_kmers);
            std::map<sgNodeID_t,uint32_t > node_match_count;
            auto rksize=read_kmers.size();
            auto longest_seed = seed_size;
            for (auto i=0;i<rksize;++i) {
                longest_seed=(longest_seed>=seed_size ? longest_seed-1:longest_seed); //decreases the length of the skipped matches in each cycle
                pmatch.node = 0;
                for (auto const &m:kmer_matches[i]) {
                    auto last_node = m.first;
                    auto last_pos = m.second;
                    auto needed_hit=i+longest_seed-k;
                    if (needed_hit>i and (needed_hit >= rksize or not in_vector(kmer_matches[needed_hit], {last_node, last_pos + needed_hit -i})))
                        continue;
                    auto last_i = i;
                    while (last_i + 1 < rksize and
                           in_vector(kmer_matches[last_i + 1], {last_node, last_pos + last_i - i + 1}))
                        ++last_i;
                    if (last_i == needed_hit and i==pmatch.read_position) {
                        pmatch.node = 0;
                    } else if (last_i > needed_hit) {
                        pmatch.node = last_node;
                        pmatch.node_position = (last_node>0 ? last_pos:last_pos+1-k);
                        pmatch.read_position = i;
                        longest_seed = pmatch.size = k + last_i - i;
                    }
                }
                if (pmatch.node != 0) {
                    private_read_perfect_matches.emplace_back(pmatch);
                }
            }
            read_perfect_matches[rid] = private_read_perfect_matches;
        }
    }
}

void LongReadsRecruiter::map(uint16_t seed_size, uint64_t first_read, uint64_t last_read) {
    if (seed_size<k) seed_size=k;
    NKmerIndex nki(sdg,k,f);
    if (last_read==0) last_read=datastore.size();
#pragma omp parallel shared(nki)
    {
        StreamKmerFactory skf(k);
        std::vector<std::pair<bool,uint64_t >> read_kmers; //TODO: type should be defined in the kmerisers
        std::vector<std::vector<std::pair<int32_t,int32_t >>> kmer_matches; //TODO: type should be defined in the index class
        PerfectMatchExtender pme(sdg,k);
        ReadSequenceBuffer sequence_reader(datastore);
        uint64_t private_mapped_reads=0;
#pragma omp for
        for (auto rid=first_read;rid<=last_read;++rid){
            std::vector<PerfectMatch> private_read_perfect_matches;
            const auto seq=std::string(sequence_reader.get_read_sequence(rid));
            read_kmers.clear();
            skf.produce_all_kmers(seq.c_str(),read_kmers);
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
                    if (pme.last_readpos-rki>=seed_size) {
                        if (!pmebp.empty()) {
                            rki = pme.last_readpos + 1 - k; //avoid extra index lookups for kmers already used once
                            pme.make_path_as_perfect_matches();
                            private_read_perfect_matches.insert(private_read_perfect_matches.end(),pme.best_path_matches.begin(),pme.best_path_matches.end());
                        }
                    }
//                    //TODO: matches shold be extended left to avoid unneeded indeterminaiton when an error occurrs in an overlap region and the new hit matches a further part of the genome.
//                    std::cout<<"rki after match extension: "<<rki<<" / "<<read_kmers.size()<<std::endl;
                }

            }
            read_perfect_matches[rid] = private_read_perfect_matches;
            if (++private_mapped_reads%5000==0) {
                sdglib::OutputLog()<<private_mapped_reads<<" reads mapped on thread "<<omp_get_thread_num()<<std::endl;
            }
        }
    }
}

//TODO: add read position and node position to match
void LongReadsRecruiter::recruit_reads(uint16_t seed_size, uint16_t seed_count, uint64_t first_read,
                                       uint64_t last_read) {
    node_reads.clear();
    node_reads.resize(sdg.nodes.size());
    if (last_read==0) last_read=datastore.size();
    for (auto rid=first_read;rid<=last_read;++rid) {
        std::map<sgNodeID_t, uint32_t> node_match_count;
        for (const auto &pmatch:read_perfect_matches[rid]) {
            if (pmatch.size>=seed_size) {
                node_match_count[pmatch.node] = node_match_count[pmatch.node] + 1;
                if (node_match_count[pmatch.node] == seed_count) {
                    node_reads[llabs(pmatch.node)].emplace_back(rid);
                }
            }
        }
    }
}
void LongReadsRecruiter::recruit_threads() {
    node_threads.clear();
    node_threads.resize(sdg.nodes.size());
    for (auto rid=1;rid<=read_threads.size();++rid)
        for (const auto &match:read_threads[rid])
                    node_threads[llabs(match.node)].emplace_back(rid);
}

std::vector<NodePosition> LongReadsRecruiter::endmatches_to_positions(uint64_t rid, int32_t end_size, uint16_t matches) {
    std::vector<NodePosition> positions;
    float MAX_DIFF=.1;
    std::map<sgNodeID_t, int> node_count;
    for (auto &pm : read_perfect_matches[rid]) node_count[pm.node]+=1;
    PerfectMatch start_first,start_second,end_first,end_second;
    for (auto &n:node_count) {
        if (n.second>=matches) {
            int64_t nsize=sdg.get_node_size(n.first);
            //std::cout<<std::endl<<"Matches vs. node "<<n.first<< " ( "<<nsize<<"bp )"<<std::endl;
            auto mi=read_perfect_matches[rid].begin();
            for (auto i=0;i<n.second;++i){//iterate thorugh the node's matches TODO: only iterate through matches in ends? a spurious middle match can be a disaster
                while (mi->node!=n.first) ++mi; // get to the i-th match to this node
                //std::cout<<mi->read_position<<" -> "<<mi->node_position;
                if (i==0) {
                    start_first=*mi;
                    //std::cout<<" SF ";
                }
                if (i==matches-1) {
                    start_second=*mi;
                    //std::cout<<" SS ";
                }
                if (i==n.second-matches) {
                    end_first=*mi;
                    //std::cout<<" EF ";
                }
                if (i==n.second-1) {
                    end_second=*mi;
                    //std::cout<<" ES ";
                }
                //std::cout<<std::endl;
                ++mi;
            }

            //std::cout<<start_first.node_position<<" "<<start_second.node_position<<" : "<<end_first.node_position<<" "<<end_second.node_position<<"  ( node: "<<nsize<<"  last_match:"<<read_perfect_matches[rid].back().read_position<<" )"<<std::endl;
            int32_t nstart=INT32_MIN,nend=INT32_MIN;
            //std::cout<<start_second.node_position<<"<="<<end_size<<"?"<<std::endl;
            if (start_second.node_position<=end_size) nstart=start_first.read_position-start_first.node_position;
            //std::cout<<end_first.node_position<<">="<<nsize-end_size<<"?"<<std::endl;
            if (end_first.node_position>=nsize-end_size) nend=end_second.read_position-end_second.node_position+nsize;
            //std::cout<<"from matches --> "<<nstart<<" - "<<nend<<std::endl;
            //skip if no end detected
            if (nstart==INT32_MIN and nend==INT32_MIN) continue;

            //skip if only start and end should be there
            if (nend==INT32_MIN){
                if (nstart+nsize<=read_perfect_matches[rid].back().read_position) continue;
                nend=nstart+nsize;
            }
            //skip if only end and start should be there
            else if (nstart==INT32_MIN){
                if (nend-nsize>0) continue;
                nstart=nend-nsize;
            }
            //skip if node size is wrong
            else if (llabs(nend-nstart-nsize)/nsize>MAX_DIFF) continue;
            //std::cout<<"adjusted to --> "<<nstart<<" - "<<nend<<std::endl;
            positions.emplace_back(n.first,nstart,nend);

        }
    }
    std::sort(positions.begin(),positions.end());
    return positions;
}

void LongReadsRecruiter::thread_reads(uint32_t end_size, uint16_t matches) {
    read_threads.clear();
    read_threads.resize(read_perfect_matches.size());
    for (auto rid=1;rid<read_perfect_matches.size();++rid){
        if (read_perfect_matches[rid].empty()) continue;
        //std::cout<<std::endl<<"analysing read "<<rid<<" ( matches from "<<read_perfect_matches[rid].front().read_position<<" to "<<read_perfect_matches[rid].back().read_position<<" )"<<std::endl;
        read_threads[rid]=endmatches_to_positions(rid,end_size,matches);
    }
}

DistanceGraph LongReadsRecruiter::dg_from_threads(bool multi_link) {
    DistanceGraph dg(sdg,"LRRthreads");
    for (auto rid=1;rid<read_threads.size();++rid) {
        auto &pos=read_threads[rid];
        if (pos.empty()) continue;
        if (multi_link) {
            for (auto i = 0; i < pos.size() - 1; ++i) {
                for (auto j = i + 1; j < pos.size(); ++j) {
                    dg.add_link(-pos[i].node, pos[j].node, pos[j].start - pos[i].end, {SupportType::LongRead, 0, static_cast<uint64_t>(rid)});
                }
            }
        }
        else {
            for (auto i = 0; i < pos.size() - 1; ++i) {
                dg.add_link(-pos[i].node, pos[i + 1].node, pos[i + 1].start - pos[i].end, {SupportType::LongRead, 0, static_cast<uint64_t>(rid)});
            }
        }
    }
    return dg;
}



//*********************** WELL INTENTIONED WARNING ****************************
// Thread And Pop is a nugget of what it could potentially become.
// It is unreliable and buggy, and so will be for a while.It works where it has
// been tested, and I am working on making it better. For now, it still is my
// pet project.
//
// bj
//*****************************************************************************

class AggregatedHits{
public:
    AggregatedHits(const std::vector<PerfectMatch> &_read_pms,uint start):read_pms(_read_pms){
        start_idx=start;
        for(end_idx=start;end_idx+1<read_pms.size();++end_idx) {
            auto & l = read_pms[end_idx];
            auto & n = read_pms[end_idx+1];
            if (n.node != l.node or n.read_position < l.read_position or l.node_position < l.node_position) break;
        }
    }

    uint hit_count() const { return end_idx-start_idx+1;}

    sgNodeID_t node() const { return read_pms[start_idx].node;}

    int32_t rstart() const {return read_pms[start_idx].read_position;}

    int32_t rend() const {return read_pms[end_idx].read_position+read_pms[end_idx].size;}

    int32_t nstart() const {return read_pms[start_idx].node_position;}

    int32_t nend() const {return read_pms[end_idx].node_position+read_pms[end_idx].size;}

    bool can_join(const AggregatedHits &other) const {
        if (node() == other.node()) {
            if (rend() < other.rstart() and nend() < other.nend()) return true;
            if (other.rend() < rstart() and other.nend() < nend()) return true;
        }
        return false;
    }

    uint start_idx;
    uint end_idx;
    const std::vector<PerfectMatch> &read_pms;
};

std::vector<AggregatedHits> aggregate_hits(const std::vector<PerfectMatch>&prm) {
    std::vector<AggregatedHits> ahs;
    for(auto next_hit=0;next_hit<prm.size();){
        ahs.emplace_back(AggregatedHits(prm,next_hit));
        next_hit=ahs.back().end_idx+1;
    }
    return ahs;
}

std::vector<sgNodeID_t> LongReadsRecruiter::path_fw(seqID_t read_id, sgNodeID_t node) const {
    //check for obvious no-path conditions:
    if (read_paths[llabs(read_id)].size()<2) return {};
    //get read paths in the right orientation - leaves r2p empty if not using pairs
    std::vector<sgNodeID_t> r1p;
    if (read_id > 0) {
        r1p=read_paths[read_id];
    }
    else {
        sdglib::reverse_path(r1p,read_paths[-read_id]); //no read pair, as the pair comes "BEFORE"
    }

    //discards up to node in r1p
    auto r1it=r1p.cbegin();
    while (r1it<r1p.cend() and *r1it!=node) ++r1it;
    if (r1it==r1p.cend()) return {};
    std::vector<sgNodeID_t> path_fw;
    path_fw.insert(path_fw.end(),++r1it,r1p.cend());

    return path_fw;
}

std::vector<std::vector<sgNodeID_t> > LongReadsRecruiter::all_paths_fw(sgNodeID_t node) const {
    std::vector<std::vector<sgNodeID_t> > r;
    std::unordered_set<seqID_t> rids;
    for (auto rid:node_paths[llabs(node)] ){
        if (node<0) rid=-rid;
        r.emplace_back(path_fw(rid,node));
        if (r.back().empty()) r.pop_back();

    }
    return r;
}

//def rudimentary_tap(rid):
//    #transform consecutive runs of ascending hits to the same node into single long-runs
//    prm=[x for x in lrr.read_perfect_matches[rid]]
//    ahs=aggregate_hits(prm)
//
//    #identify runs of just one hit, evaluate if they represent a transition or a switch
//    switches=[]
//    if ahs[0].size()==1: switches.append(0)
//    if ahs[-1].size()==1: switches.append(len(prm)-1)
//    for i in range(1,len(ahs)-1):
//        if ahs[i].size()==1 and ahs[i-1].can_join(ahs[i+1]):
//            #print(ahs[i]," identified as a 1-hit switch")
//            for x in range(ahs[i].start_index,ahs[i].end_index+1): switches.append(x)
//    prm=[prm[x] for x in range(len(prm)) if x not in switches]
//    ahs=aggregate_hits(prm)
//    return ahs

std::vector<sgNodeID_t> rudimentary_tap(const std::vector<PerfectMatch>&prm){
    std::vector<PerfectMatch> filtered_prm;
    auto ahs=aggregate_hits(prm);

    for (auto i=0;i<ahs.size();++i){
        if (i==0 or i==ahs.size()-1){
            if (ahs[i].hit_count()==1) continue; //Discards starting and ending single-hits, may be valid but not trustworthy
        }
        else {
            if (ahs[i].hit_count()==1 and ahs[i-1].can_join(ahs[i+1])) continue; //Discards 1-hit switches
        }
        for (auto mi=ahs[i].start_idx;mi<=ahs[i].end_idx;++mi) filtered_prm.emplace_back(prm[mi]);
    }
    std::vector<sgNodeID_t> p;
    for (auto &ah: aggregate_hits(filtered_prm)) p.emplace_back(ah.node());
    return p;
}

void LongReadsRecruiter::thread_and_pop() {
    read_paths.clear();
    read_paths.resize(read_perfect_matches.size());
    node_paths.clear();
    node_paths.resize(sdg.nodes.size());
    std::vector<uint32_t> pcount(sdg.nodes.size());
    for (auto rid=1;rid<read_perfect_matches.size();++rid) {
        read_paths[rid]=rudimentary_tap(read_perfect_matches[rid]);
        for (auto const & n:read_paths[rid]) ++pcount[abs(n)];
    }
    for(auto i=1;i<pcount.size();++i) node_paths[i].reserve(pcount[i]);
    for (auto rid=1;rid<read_perfect_matches.size();++rid){
        for (auto const & n:read_paths[rid]) node_paths[abs(n)].emplace_back( n>0 ? rid : -rid);
    }
}
