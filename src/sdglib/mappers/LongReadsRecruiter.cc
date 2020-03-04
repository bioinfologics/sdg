//
// Created by Bernardo Clavijo (EI) on 16/01/2020.
//

#include "LongReadsRecruiter.hpp"
#include <sdglib/indexers/NKmerIndex.hpp>

LongReadsRecruiter::LongReadsRecruiter(const SequenceDistanceGraph &_sdg, const LongReadsDatastore &datastore,uint8_t _k=25, uint16_t _f=50):
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

DistanceGraph LongReadsRecruiter::dg_from_threads() {
    DistanceGraph dg(sdg);
    for (auto rid=1;rid<read_threads.size();++rid) {
        auto &pos=read_threads[rid];
        if (pos.empty()) continue;
        for (auto i = 0; i < pos.size() - 1; ++i) {
            //std::cout<<"dg.add_link("<<-pos[i].node<<","<<pos[i+1].node<<","<<pos[i+1].start-pos[i].end<<",{SupportType::LongRead,0,"<<static_cast<uint64_t>(rid)<<"})"<<std::endl;
            dg.add_link(-pos[i].node, pos[i + 1].node, pos[i + 1].start - pos[i].end,
                        {SupportType::LongRead, 0, static_cast<uint64_t>(rid)});
        }
    }
}