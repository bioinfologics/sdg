//
// Created by Luis Yanes (EI) on 22/03/2018.
//

#include <cmath>
#include <sglib/mappers/LongReadMapper.hpp>
#include <sglib/logger/OutputLog.hpp>

#include <sglib/utilities/omp_safe.hpp>
#include <atomic>
#include <sglib/utilities/io_helpers.hpp>

const bsgVersion_t LongReadMapper::min_compat = 0x0001;

void LongReadHaplotypeMappingsFilter::set_read(uint64_t _read_id) {
    read_id=_read_id;
    mappings=lorm.get_raw_mappings_from_read(read_id);
}

void LongReadHaplotypeMappingsFilter::generate_haplotypes_from_linkedreads(float min_tn) {

}

std::array<uint64_t,3> LongReadHaplotypeMappingsFilter::haplotype_coverage_dist(
        const std::vector<sgNodeID_t> &haplotype) {
    return {0,0,0};
}

void LongReadHaplotypeMappingsFilter::score_coverage(float weight) {

}

void LongReadHaplotypeMappingsFilter::score_window_winners(int win_size, int win_step, float weight) {

}

void LongReadHaplotypeMappingsFilter::rank(uint64_t read_id, float coverage_weight, float winners_weight) {

}

std::vector<LongReadMapping> LongReadMapper::get_raw_mappings_from_read(uint64_t read_id) const {
    std::vector<LongReadMapping> r;
    if (first_mapping[read_id]!=-1) {
        for (auto i = first_mapping[read_id];mappings[i].read_id==read_id;++i) r.emplace_back(mappings[i]);
    }
    return std::move(r);
}

void LongReadMapper::get_all_kmer_matches(std::vector<std::vector<std::pair<int32_t, int32_t>>> & matches, std::vector<std::pair<bool, uint64_t>> & read_kmers){
    if (matches.size() < read_kmers.size()) matches.resize(read_kmers.size());
    uint64_t no_match=0,single_match=0,multi_match=0; //DEBUG
    for (auto i=0;i<read_kmers.size();++i){
        matches[i].clear();
        auto first = assembly_kmers.find(read_kmers[i].second);
        for (auto it = first; it != assembly_kmers.end() && it->kmer == read_kmers[i].second; ++it) {
            int32_t offset=it->offset; //so far, this is +1 and the sign indicate direction of kmer in contig
            sgNodeID_t node=it->contigID; //so far, this is always positive
            if (read_kmers[i].first != (offset>0) ) {
                node=-node;
                offset=( (int) sg.nodes[std::llabs(it->contigID)].sequence.size() ) - std::abs(offset);
            }
            else offset=std::abs(offset)-1;
            matches[i].emplace_back(node, offset);
        }
        if (matches[i].empty()) ++no_match; //DEBUG
        else if (matches[i].size()==1) ++single_match; //DEBUG
        else ++multi_match; //DEBUG
    }
//    std::cout<<"From get_all_kmer_matches with "<< read_kmers.size() <<": "<<no_match<<" ("<< (no_match*100/read_kmers.size()) << "%) no match   "
//            <<single_match<<" ("<< (single_match*100/read_kmers.size()) << "%) single match   "
//            <<multi_match<<" ("<< (multi_match*100/read_kmers.size()) << "%) multi match"<<std::endl;
}

std::set<sgNodeID_t> LongReadMapper::window_candidates(std::vector<std::vector<std::pair<int32_t, int32_t>>> & matches, uint32_t read_kmers_size){
    //beware you need the size! otherwise there is a size from the matches and not from the reads!
    std::set<sgNodeID_t> candidates;
    std::map<int32_t,uint32_t> candidate_votes;
    std::map<int32_t,uint32_t> candidate_multivotes;
    for (auto i=0;i<read_kmers_size;++i) {
        for (auto m:matches[i]) ++candidate_votes[m.first];
        if (matches[i].size()>1) for (auto m:matches[i]) ++candidate_multivotes[m.first];
    }
    //std::cout<<"From window_candidates, hit nodes:"<<std::endl;
    //for (auto cv:candidate_votes) std::cout<<"   "<<cv.first<<": "<<cv.second<<" ("<<candidate_multivotes[cv.first]<<" multi "<<cv.second-candidate_multivotes[cv.first]<<" single)"<<std::endl;
    //std::cout<<std::endl;
    for (auto cv:candidate_votes) if (cv.second>50) candidates.insert(cv.first);

    return candidates;
}

std::vector<LongReadMapping> LongReadMapper::alignment_blocks(uint32_t readID, std::vector<std::vector<std::pair<int32_t, int32_t>>> & matches,
                                                              uint32_t read_kmers_size,
                                                              std::set<sgNodeID_t> &candidates) {
//    std::cout<<"Alignment block search starting with candidates:";
//    for (auto c:candidates) std::cout<<" "<<c;
//    std::cout<<std::endl;
    std::vector<LongReadMapping> blocks;
    //transform to matches per candidate;
    std::map<int32_t , std::vector<std::pair<int32_t,int32_t>>> candidate_hits;

    for (auto p=0;p<read_kmers_size;++p){
        for (auto m:matches[p]) {
            if (candidates.count(m.first)){
                //std::cout<<"adding hit for candidate "<<m.first<<" "<<p<<" -> "<<m.second<<std::endl;
                candidate_hits[m.first].emplace_back(p,m.second);
            }
        }
    }

    for (auto &ch:candidate_hits){
        int start_p=0;
        int start_t=0;
        int last_delta=0;
        int last_t=-1000;
        int last_p=-1000;
        int chain=-1;
        auto &chits=ch.second;
        std::vector<bool> used(ch.second.size());
        //std::sort(chits.begin(),chits.end());
        //Filter hits on overused positions
        for (auto starti=0 ; starti<chits.size() ; ++starti){
            auto endi=starti;
            while (endi<chits.size()-1 and chits[starti].first==chits[endi+1].first) ++endi;
            bool discard=false;
            if (endi>starti+4) discard=true; //more than 5 hits to candidate ->discard
            else {
                for (auto x=starti;x<endi;++x) {
                    if (chits[x+1].second-chits[x].second<max_jump) discard=true;
                }
            }
            if (discard) for (auto x=starti;x<=endi;++x) used[x]=true;
        }

        for (auto starti=0 ; starti<chits.size() ; ++starti){
            if (used[starti]) continue;
            chain=1;
            start_p=chits[starti].first;
            start_t=chits[starti].second;
            last_p=chits[starti].first;
            last_t=chits[starti].second;
            last_delta=last_t-last_p;
            used[starti]=true;
            //std::cout<<"New chain to "<<ch.first<<" started with "<<start_p<<"->"<<start_t<<std::endl;
            for (auto i=starti+1;i<chits.size() and chits[i].first-last_p<=max_jump;++i){
                auto &h=chits[i];
                //if not in chain, continue;
                if (used[i] or last_p > h.first or last_t > h.second or h.second-last_t > max_jump
                    or llabs(last_delta-(h.second-h.first))>max_delta_change) {
                    //std::cout<<" skipping "<<h.first<<"->"<<h.second<<std::endl;
                    continue;
                }
                ++chain;
                last_p=h.first;
                last_t=h.second;
                last_delta=last_t-last_p;
                used[i]=true;
                //if (start_p==61 and start_t==3621) std::cout<<" chain is now length "<<chain<<" after adding "<<h.first<<"->"<<h.second<<" (delta "<<h.first-h.second<<")"<<std::endl;
            }
            if (chain>=min_chain and last_p-start_p>=min_size) {
                LongReadMapping b;
                b.read_id=readID;
                b.qStart=start_p;
                b.qEnd=last_p;
                b.node=ch.first;
                b.nStart=start_t;
                b.nEnd=last_t;
                b.score=chain;
                blocks.push_back(b);
            }

        }
    }
    std::sort(blocks.begin(),blocks.end());
    return blocks;
}

std::vector<LongReadMapping> LongReadMapper::filter_blocks(std::vector<LongReadMapping> &blocks,
                                                           std::vector<std::vector<std::pair<int32_t, int32_t>>> &matches,
                                                           uint32_t read_kmers_size) {
    std::vector<LongReadMapping> fblocks;
    for (auto &b:blocks){
        bool looser=false;
        for (auto &ob:blocks){
            if (ob==b) continue;
            if (ob.qStart<=b.qStart and ob.qEnd>=b.qEnd and ob.score>b.score) {
                looser=true;
                break;
            }
        }
        if (!looser) fblocks.emplace_back(b);
    }
    return fblocks;
}

void LongReadMapper::map_reads(std::unordered_set<uint32_t> readIDs, std::string detailed_log) {
    update_graph_index();
    std::vector<std::vector<LongReadMapping>> thread_mappings(omp_get_max_threads());
    std::atomic<uint32_t > num_reads_done(0);
    std::atomic<uint64_t > window_low_score(0),window_close_second(0),window_hit(0);
    std::atomic<uint64_t > no_matches(0),single_matches(0),multi_matches(0);
    std::ofstream dl;
    if (!detailed_log.empty()) dl.open(detailed_log);
#pragma omp parallel
    {
        StreamKmerFactory skf(k);
        std::vector<std::pair<bool, uint64_t>> read_kmers;

        auto & private_results=thread_mappings[omp_get_thread_num()];
        BufferedSequenceGetter sequenceGetter(datastore);
        std::vector<std::vector<std::pair<int32_t, int32_t>>> node_matches; //node, offset
        const char * query_sequence_ptr;
#pragma omp for
        for (uint32_t readID = 1; readID < datastore.size(); ++readID) {

            if (readID%1000 == 0) {
                num_reads_done+=1000;
#pragma omp critical (lrmap_progress)
                {
                    sglib::OutputLog() << num_reads_done << " / " << datastore.size() << " reads mapped" << std::endl;
                }
            }

            if ((!readIDs.empty() and readIDs.count(readID)==0 /*Read not in selected set*/)
                or
                (readIDs.empty()/*No selection*/)) {
                continue;
            }

            //========== 1. Get read sequence, kmerise, get all matches ==========
            query_sequence_ptr = sequenceGetter.get_read_sequence(readID);
            read_kmers.clear();
            skf.produce_all_kmers(query_sequence_ptr, read_kmers);

            if ( read_kmers.size()< 2 * min_size) {
                continue;
            }

            get_all_kmer_matches(node_matches,read_kmers);

            //========== 2. Find match candidates in fixed windows ==========

            auto candidates = window_candidates(node_matches,read_kmers.size());

            //========== 3. Create alignment blocks from candidates ==========

            auto blocks = alignment_blocks(readID,node_matches,read_kmers.size(),candidates);

            //========== 4. Construct mapping path ==========
            if (blocks.empty()) ++no_matches;
            else if (blocks.size()==1) ++single_matches;
            else ++multi_matches;
            //TODO: align blocks that occupy the same space as longer/better blocks should be thrown away.

            //auto fblocks=filter_blocks(blocks,node_matches,read_kmers.size());
            auto fblocks=blocks;

//            std::cout<<"After FILTERING:"<<std::endl;
//            for (auto b:fblocks)
//                std::cout << "Target: " << b.node << " (" << sg.nodes[llabs(b.node)].sequence.size() << " bp)  "
//                          << b.qStart << ":" << b.qEnd << " -> " << b.nStart << ":" << b.nEnd
//                          << " (" << b.score << " chained hits, " << b.qEnd - b.qStart + k << "bp, "
//                          << b.score * 100 / (b.qEnd - b.qStart) << "%)"
//                          << std::endl;

            private_results.insert(private_results.end(),fblocks.begin(),fblocks.end());


        }
    }
    //TODO: report read and win coverage by winners vs. by loosers
//    sglib::OutputLog()<<"Read window results:    "<<window_low_score<<" low score    "<<window_close_second<<" close second    "<<window_hit<<" hits"<<std::endl;
    sglib::OutputLog()<<"Read results:    "<<no_matches<<" no match    "<<single_matches<<" single match    "<<multi_matches<<" multi matches"<<std::endl;
    for (auto & threadm : thread_mappings) {
        mappings.insert(mappings.end(),threadm.begin(),threadm.end());
    }
    sglib::OutputLog() << "Updating mapping indexes" <<std::endl;
    update_indexes();
}

void LongReadMapper::update_indexes() {
    std::sort(mappings.begin(),mappings.end());
    first_mapping.clear();
    first_mapping.resize(datastore.size(),-1);
    for (int64_t i=0;i<mappings.size();++i){
        auto &m=mappings[i];
        if (first_mapping.size()<=m.read_id) first_mapping.resize(m.read_id+1,-1);
        if (first_mapping[m.read_id]==-1) first_mapping[m.read_id]=i;
    }
    reads_in_node.clear();
    reads_in_node.resize(sg.nodes.size());
    for (auto &rm:filtered_read_mappings) {
        for (auto &m:rm) {
            if (reads_in_node[std::abs(m.node)].empty() or reads_in_node[std::abs(m.node)].back()!=m.read_id)
                reads_in_node[std::abs(m.node)].push_back(m.read_id);
        }
    }
    sglib::OutputLog() << "Finished with " << filtered_read_mappings.size() << " mappings" << std::endl;
}

LongReadMapper::LongReadMapper(SequenceGraph &sg, LongReadsDatastore &ds, uint8_t k)
        : sg(sg), k(k), datastore(ds), assembly_kmers(k) {

    reads_in_node.resize(sg.nodes.size());
}

LongReadMapper::~LongReadMapper() {}

void LongReadMapper::read(std::string filename) {
    std::ifstream inf(filename, std::ios_base::binary);
    read(inf);
}

void LongReadMapper::read(std::ifstream &inf) {
    sglib::OutputLog() << "Reading long read mappings" << std::endl;
    auto mapSize(mappings.size());
    inf.read(reinterpret_cast<char *>(&k), sizeof(k));
    inf.read(reinterpret_cast<char *>(&mapSize), sizeof(mapSize));
    mappings.resize(mapSize);
    inf.read(reinterpret_cast<char*>(mappings.data()), mappings.size()*sizeof(LongReadMapping));
    sglib::OutputLog() << "Updating read mapping indexes!" << std::endl;
    update_indexes();
    sglib::OutputLog() << "Done!" << std::endl;
}

void LongReadMapper::write(std::string filename) {
    std::ofstream ofs(filename, std::ios_base::binary);
    write(ofs);
}

void LongReadMapper::write(std::ofstream &ofs) {
    sglib::OutputLog() << "Dumping long read mappings" << std::endl;
    auto mapSize(mappings.size());
    ofs.write(reinterpret_cast<const char *>(&k), sizeof(k));
    ofs.write(reinterpret_cast<const char *>(&mapSize), sizeof(mapSize));
    ofs.write(reinterpret_cast<const char*>(mappings.data()), mappings.size()*sizeof(LongReadMapping));
    sglib::OutputLog() << "Done!" << std::endl;
}

void LongReadMapper::write_filtered_mappings(std::string filename) {
    std::ofstream ofs(filename, std::ios_base::binary);
    sglib::write_flat_vectorvector(ofs,filtered_read_mappings);
}

void LongReadMapper::read_filtered_mappings(std::string filename) {
    std::ifstream ifs(filename, std::ios_base::binary);
    sglib::read_flat_vectorvector(ifs,filtered_read_mappings);
}

LongReadMapper LongReadMapper::operator=(const LongReadMapper &other) {
    if (&sg != &other.sg and &datastore != &other.datastore) { throw ("Can only copy paths from the same SequenceGraph"); }
    if (&other == this) {
        return *this;
    }
    mappings = other.mappings;
    update_indexes();
    return *this;
}

void LongReadMapper::update_graph_index() {
    std::cout<<"updating index with k="<<std::to_string(k)<<std::endl;
    assembly_kmers=NKmerIndex(k);
    assembly_kmers.generate_index(sg);
}

void LongReadMapper::filter_mappings_with_linked_reads(const LinkedReadMapper &lrm, uint32_t min_size, float min_tnscore) {
    if (lrm.tag_neighbours.empty()) {
        sglib::OutputLog()<<"Can't filter mappings because there are no tag_neighbours on the LinkedReadMapper"<<std::endl;
        return;
    }
    sglib::OutputLog()<<"Filtering mappings"<<std::endl;
    filtered_read_mappings.clear();
    filtered_read_mappings.resize(datastore.size());
    BufferedSequenceGetter lrbsg(datastore);
    //create a collection of mappings for every read.
    uint64_t last_mapindex=0;
    uint64_t rshort=0,runmapped=0,rnocovset=0,rprocessed=0,rnoset=0;
    for (uint64_t rid=0;rid<datastore.size();++rid) {
        while (last_mapindex<mappings.size() and mappings[last_mapindex].read_id<rid) ++last_mapindex;
        auto r=filter_mappings_with_linked_reads(lrm,lrbsg,min_size,min_tnscore,rid,last_mapindex);
        switch(r){
            case MappingFilterResult::Success:
                ++rprocessed;
                break;
            case MappingFilterResult::TooShort:
                ++rshort;
                break;
            case MappingFilterResult::NoMappings:
                ++runmapped;
                break;
            case MappingFilterResult::NoReadSets:
                ++rnoset;
                break;
            case MappingFilterResult::LowCoverage:
                ++rnocovset;
                break;
        }
    }
    sglib::OutputLog()<<"too short: "<<rshort<<"  no mappings: "<<runmapped<<"  no sets: "<<rnoset<<"  no cov1>50%: "<<rnocovset<<" processed: "<<rprocessed<<std::endl;
}

MappingFilterResult LongReadMapper::filter_mappings_with_linked_reads(const LinkedReadMapper &lrm, BufferedSequenceGetter &lrbsg, uint32_t min_size, float min_tnscore, uint64_t readID, uint64_t offset_hint) {
    uint64_t last_mapindex=offset_hint;
    const char* seq_ptr;
    seq_ptr = lrbsg.get_read_sequence(readID);
    std::string seq(seq_ptr);
    if (seq.size()<min_size) return MappingFilterResult::TooShort;
    auto rsize=seq.size();
    while (last_mapindex>0 and mappings[last_mapindex].read_id>readID) --last_mapindex;
    while (last_mapindex<mappings.size() and mappings[last_mapindex].read_id<readID) ++last_mapindex;
    if (last_mapindex>=mappings.size() or mappings[last_mapindex].read_id!=readID) return MappingFilterResult::NoMappings;
    std::vector<LongReadMapping> read_mappings;
    std::unordered_set<sgNodeID_t> all_nodes;
    while (last_mapindex<mappings.size() and mappings[last_mapindex].read_id==readID) {
        read_mappings.emplace_back(mappings[last_mapindex]);
        all_nodes.insert(llabs(read_mappings.back().node));
        ++last_mapindex;
    }
    //evaluate a read's mappings:
    //TODO: 1) remove middle-mappings (later)
    //2) create the neighbour sets, uniq to not evaluate same set many times over.
    std::map<std::set<sgNodeID_t>,int> all_nsets;
    for (auto n:all_nodes){
        std::set<sgNodeID_t> nset;
        for (auto m:all_nodes){
            for (auto &score:lrm.tag_neighbours[n]){
                if (score.node==m and score.score>min_tnscore) {
                    nset.insert(m);
                    break;
                }
            }
        }
        ++all_nsets[nset];
    }
    if (all_nsets.empty()) return MappingFilterResult::NoReadSets;
    std::vector<std::pair<uint64_t,std::set<sgNodeID_t>>> cov1set;

    //3) remove all sets not passing a "sum of all mappings >50% of the read"
    //   compute the 1-cov total bases for each set, if it is not 50% of the read discard too
    for (auto &nsv:all_nsets){
        std::vector<uint8_t> coverage(seq.size());
        auto &nodeset=nsv.first;
        uint64_t total_map_bp=0;
        for (auto m:read_mappings) if (nodeset.count(llabs(m.node))) total_map_bp+=m.qEnd-m.qStart;
        if (total_map_bp*2<seq.size()) continue;
        //this can be done faster by saving all starts and ends and keeping a rolling value;
        for (auto m:read_mappings) {
            if (nodeset.count(llabs(m.node))){
                for (auto i=m.qStart;i<=m.qEnd and i<rsize;++i) ++coverage[i];
            }
        }
        uint64_t cov1=std::count(coverage.begin(),coverage.end(),1);
        if (cov1*2>=seq.size())
            cov1set.emplace_back(cov1, nodeset);
    }
    if (cov1set.empty()) return MappingFilterResult::LowCoverage;
    //find set with highest 1-cov
    std::sort(cov1set.begin(),cov1set.end());
    auto &winset=cov1set.back().second;
    //TODO: 4) try to find "turn-arounds" for CCS-style reads (later)
    //TODO: 7) copy vector to read_to_mappings
    filtered_read_mappings[readID].clear();
    for (auto &m:read_mappings){
        if (winset.count(llabs(m.node))) filtered_read_mappings[readID].push_back(m);
    }
    return MappingFilterResult::Success;
}