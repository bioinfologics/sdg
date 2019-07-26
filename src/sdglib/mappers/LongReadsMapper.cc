//
// Created by Luis Yanes (EI) on 22/03/2018.
//

#include "LongReadsMapper.hpp"
#include <atomic>
#include <cmath>
#include <iomanip>      // std::setprecision
#include <iterator>
#include <sstream>
#include <sdglib/utilities/OutputLog.hpp>
#include <sdglib/utilities/omp_safe.hpp>
#include <sdglib/utilities/io_helpers.hpp>
#include <sdglib/utilities/most_common_helper.hpp>
#include <sdglib/workspace/WorkSpace.hpp>


const sdgVersion_t LongReadsMapper::min_compat = 0x0003;

std::string LongReadsMapper::ls(int level,bool recursive) {
    std::stringstream ss;
    std::string spacer(2 * level, ' ');
    ss << spacer << "Long Reads Mapper: "<<mappings.size()<<" raw mappings" << std::endl;
    return ss.str();
}

void LongReadsMapper::print_status() const {
    sdglib::OutputLog()<<"Long read mappings: "<<mappings.size()<<" mappings."<<std::endl;
}

std::vector<LongReadMapping> LongReadsMapper::get_raw_mappings_from_read(uint64_t read_id) const {
    std::vector<LongReadMapping> r;
    if (first_mapping[read_id]!=-1) {
        for (auto i = first_mapping[read_id];i<mappings.size() and mappings[i].read_id==read_id;++i) r.emplace_back(mappings[i]);
    }
    return std::move(r);
}

void LongReadsMapper::get_sat_kmer_matches(const SatKmerIndex & sat_assembly_kmers, std::vector<std::vector<std::pair<int32_t, int32_t>>> &matches, std::vector<std::pair<bool, uint64_t>> &read_kmers) {
    if (matches.size() < read_kmers.size()) matches.resize(read_kmers.size());
    uint64_t no_match=0,single_match=0,multi_match=0; //DEBUG
    for (auto i=0;i<read_kmers.size();++i){
        matches[i].clear();
        auto start(sat_assembly_kmers.beginCO(read_kmers[i].second));
        auto end(sat_assembly_kmers.endCO(read_kmers[i].second));
        for (auto it = start; it < end; ++it) {
            int32_t offset=sat_assembly_kmers.contig_offsets[it].offset; //so far, this is +1 and the sign indicate direction of kmer in contig
            sgNodeID_t node= sat_assembly_kmers.contig_offsets[it].contigID; //so far, this is always positive
            if (read_kmers[i].first != (offset>0) ) {
                node=-node;
                offset=( (int) sg.nodes[std::llabs(sat_assembly_kmers.contig_offsets[it].contigID)].sequence.size() ) - std::abs(offset);
            }
            else offset=std::abs(offset)-1;
            matches[i].emplace_back(node, offset);
        }
        if (matches[i].empty()) ++no_match; //DEBUG
        else if (matches[i].size()==1) ++single_match; //DEBUG
        else ++multi_match; //DEBUG
    }
}

void LongReadsMapper::get_all_kmer_matches(const NKmerIndex & assembly_kmers, std::vector<std::vector<std::pair<int32_t, int32_t>>> & matches, std::vector<std::pair<bool, uint64_t>> & read_kmers) {
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

void LongReadsMapper::count_candidates(std::vector<unsigned char> &candidate_counts,
                                      std::vector<std::vector<std::pair<int32_t, int32_t>>> &matches,
                                  uint32_t read_kmers_size){
    //beware you need the size! otherwise there is a size from the matches and not from the reads!
    for (auto i=0;i<read_kmers_size;++i) {
        for (auto m:matches[i]) {
            uint32_t linear_node = (m.first>0) ? m.first : sg.nodes.size()+std::abs(m.first);
            uint16_t count_so_far=candidate_counts[ linear_node ]+1;
            candidate_counts[ linear_node ] = (count_so_far > UINT8_MAX) ? UINT8_MAX : count_so_far;
        }
    }
}

std::vector<LongReadMapping> LongReadsMapper::alignment_blocks(uint32_t readID,
                                                              std::vector<std::vector<std::pair<int32_t, int32_t>>> &matches,
                                                              uint32_t read_kmers_size,
                                                              const std::vector<unsigned char> &candidate_counts) {
//    std::cout<<"Alignment block search starting with candidates:";
//    for (auto c:candidates) std::cout<<" "<<c;
//    std::cout<<std::endl;
    std::vector<LongReadMapping> blocks;
    //transform to matches per candidate;
    std::unordered_map<int32_t , std::vector<std::pair<int32_t,int32_t>>> candidate_hits;
    candidate_hits.reserve(matches.size()/2);
    for (auto p=0;p<read_kmers_size;++p){
        for (auto m:matches[p]) {
            uint32_t linear_node = (m.first>0) ? m.first : sg.nodes.size()+std::abs(m.first);
            if (candidate_counts[linear_node] > 50) {
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

std::vector<LongReadMapping> LongReadsMapper::filter_blocks(std::vector<LongReadMapping> &blocks,
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

void LongReadsMapper::map_reads(const std::unordered_set<uint32_t> &readIDs) {
    NKmerIndex nkindex;
    SatKmerIndex skindex;
    if (sat_kmer_index) skindex=SatKmerIndex(sg,k,max_index_freq);
    else nkindex=NKmerIndex(sg,k,max_index_freq);
    std::vector<std::vector<LongReadMapping>> thread_mappings(omp_get_max_threads());
    uint32_t num_reads_done(0);
    uint64_t no_matches(0),single_matches(0),multi_matches(0);
#pragma omp parallel
    {
        StreamKmerFactory skf(k);
        std::vector<std::pair<bool, uint64_t>> read_kmers;

        auto & private_results=thread_mappings[omp_get_thread_num()];
        ReadSequenceBuffer sequenceGetter(datastore);
        std::vector<std::vector<std::pair<int32_t, int32_t>>> node_matches; //node, offset
        const char * query_sequence_ptr;

        std::vector<unsigned char> candidate_counts(sg.nodes.size()*2);

#pragma omp for schedule(static,1000) reduction(+:no_matches,single_matches,multi_matches,num_reads_done)
        for (uint32_t readID = 1; readID < datastore.size(); ++readID) {
            if (++num_reads_done%1000 == 0) {
                sdglib::OutputLog() << "Thread #"<<omp_get_thread_num() <<" processing its read #" << num_reads_done << std::endl;
            }

            if ((!readIDs.empty() and readIDs.count(readID)==0 /*Read not in selected set*/)) {
                continue;
            }
            if (datastore.read_to_fileRecord[readID].record_size< 2 * min_size ) continue;
            //========== 1. Get read sequence, kmerise, get all matches ==========
            query_sequence_ptr = sequenceGetter.get_read_sequence(readID);
            read_kmers.clear();
            skf.produce_all_kmers(query_sequence_ptr, read_kmers);

            if ( read_kmers.size()< 2 * min_size) {
                continue;
            }

            if (sat_kmer_index) {
                get_sat_kmer_matches(skindex,node_matches, read_kmers);
            } else {
                get_all_kmer_matches(nkindex,node_matches,read_kmers);
            }

            //========== 2. Find match candidates in fixed windows ==========

            count_candidates(candidate_counts, node_matches,read_kmers.size());

            //========== 3. Create alignment blocks from candidates ==========

            auto blocks(alignment_blocks(readID,node_matches,read_kmers.size(), candidate_counts));
            std::memset(candidate_counts.data(), 0, candidate_counts.size());

            //========== 4. Construct mapping path ==========
            if (blocks.empty()) ++no_matches;
            else if (blocks.size()==1) ++single_matches;
            else ++multi_matches;
            //TODO: align blocks that occupy the same space as longer/better blocks should be thrown away.

            //auto fblocks=filter_blocks(blocks,node_matches,read_kmers.size());
            const auto &fblocks = blocks;

//            std::cout<<"READ " << readID << " mappings, after FILTERING:"<<std::endl;
//            for (auto b:fblocks)
//                std::cout << "READ\tTarget: " << b.node << " (" << sdg.nodes[llabs(b.node)].sequence.size() << " bp)  "
//                          << b.qStart << ":" << b.qEnd << " -> " << b.nStart << ":" << b.nEnd
//                          << " (" << b.score << " chained hits, " << b.qEnd - b.qStart + k << "bp, "
//                          << b.score * 100 / (b.qEnd - b.qStart) << "%)"
//                          << std::endl;

            private_results.insert(private_results.end(),fblocks.begin(),fblocks.end());
        }
    }
    //TODO: report read and win coverage by winners vs. by loosers
//    sdglib::OutputLog()<<"Read window results:    "<<window_low_score<<" low score    "<<window_close_second<<" close second    "<<window_hit<<" hits"<<std::endl;
    sdglib::OutputLog()<<"Read results:    "<<no_matches<<" no match    "<<single_matches<<" single match    "<<multi_matches<<" multi matches"<<std::endl;
    for (auto & threadm : thread_mappings) {
        mappings.insert(mappings.end(),threadm.begin(),threadm.end());
    }
    sdglib::OutputLog() << "Updating mapping indexes" <<std::endl;
    update_indexes();
}

std::vector<LongReadMapping> LongReadsMapper::map_sequence(const NKmerIndex &nkindex, const char * query_sequence_ptr, sgNodeID_t seq_id) {
    std::vector<LongReadMapping> private_results;

    StreamKmerFactory skf(k);
    std::vector<std::pair<bool, uint64_t>> read_kmers;
    std::vector<std::vector<std::pair<int32_t, int32_t>>> node_matches; //node, offset
    read_kmers.clear();

    std::vector<unsigned char> candidate_counts(sg.nodes.size()*2);
    skf.produce_all_kmers(query_sequence_ptr, read_kmers);
    get_all_kmer_matches(nkindex,node_matches,read_kmers);

    //========== 2. Find match candidates in fixed windows ==========

    count_candidates(candidate_counts, node_matches,read_kmers.size());

    //========== 3. Create alignment blocks from candidates ==========

    return alignment_blocks(seq_id,node_matches,read_kmers.size(), candidate_counts);

}

void LongReadsMapper::update_indexes() {
    if (!std::is_sorted(mappings.begin(), mappings.end())) {
        sdglib::OutputLog() << "Sorting mappings" << std::endl;
        sdglib::sort(mappings.begin(), mappings.end());
    }
    reads_in_node.clear();
    reads_in_node.resize(sg.nodes.size());
    first_mapping.clear();
    first_mapping.resize(datastore.size(),-1);
    for (int64_t i=0;i<mappings.size();++i){
        auto &m=mappings[i];
        if (first_mapping.size()<=m.read_id) first_mapping.resize(m.read_id+1,-1);
        if (first_mapping[m.read_id]==-1) first_mapping[m.read_id]=i;
        if (reads_in_node[std::abs(m.node)].empty() or reads_in_node[std::abs(m.node)].back()!=m.read_id)
            reads_in_node[std::abs(m.node)].push_back(m.read_id);
    }
}

LongReadsMapper::LongReadsMapper(const SequenceDistanceGraph &_sdg, const LongReadsDatastore &ds, uint8_t k, bool sat_index)
        : sg(_sdg), k(k), datastore(ds), sat_kmer_index(sat_index) {
    reads_in_node.resize(sg.nodes.size());
}

void LongReadsMapper::read(std::string filename) {
    std::ifstream inf(filename, std::ios_base::binary);
    read(inf);
}

void LongReadsMapper::read(std::ifstream &inf) {
    sdglib::OutputLog() << "Reading long read mappings" << std::endl;

    inf.read(reinterpret_cast<char *>(&k), sizeof(k));
    inf.read(reinterpret_cast<char *>(&max_index_freq), sizeof(max_index_freq));
    inf.read(reinterpret_cast<char *>(&min_size), sizeof(min_size));
    inf.read(reinterpret_cast<char *>(&min_chain), sizeof(min_chain));
    inf.read(reinterpret_cast<char *>(&max_jump), sizeof(max_jump));
    inf.read(reinterpret_cast<char *>(&max_delta_change), sizeof(max_delta_change));

    sdglib::read_flat_vector(inf, mappings);

    sdglib::OutputLog() << "Updating read mapping indexes!" << std::endl;
    update_indexes();
    sdglib::OutputLog() << "Done!" << std::endl;
}

void LongReadsMapper::write(std::string filename) {
    std::ofstream ofs(filename, std::ios_base::binary);
    write(ofs);
}

void LongReadsMapper::write(std::ofstream &ofs) {
    sdglib::OutputLog() << "Dumping long read mappings" << std::endl;

    ofs.write(reinterpret_cast<const char *>(&k), sizeof(k));
    ofs.write(reinterpret_cast<const char *>(&max_index_freq), sizeof(max_index_freq));
    ofs.write(reinterpret_cast<const char *>(&min_size), sizeof(min_size));
    ofs.write(reinterpret_cast<const char *>(&min_chain), sizeof(min_chain));
    ofs.write(reinterpret_cast<const char *>(&max_jump), sizeof(max_jump));
    ofs.write(reinterpret_cast<const char *>(&max_delta_change), sizeof(max_delta_change));

    sdglib::write_flat_vector(ofs, mappings);
    sdglib::OutputLog() << "Done!" << std::endl;
}

void LongReadsMapper::write_read_paths(std::string filename) {
    std::ofstream ofs(filename, std::ios_base::binary);
    sdglib::write_flat_vectorvector(ofs, read_paths);
}

void LongReadsMapper::read_read_paths(std::string filename) {
    std::ifstream ifs(filename, std::ios_base::binary);
    sdglib::read_flat_vectorvector(ifs,read_paths);
}

LongReadsMapper& LongReadsMapper::operator=(const LongReadsMapper &other) {
    if (&other == this) {
        return *this;
    }
    if (&sg != &other.sg and &datastore != &other.datastore) { throw std::runtime_error("Can only LongReadsMappers from the same SequenceDistanceGraph and LongReadsDatastore"); }
    sat_kmer_index = other.sat_kmer_index;
    k = other.k;
    mappings = other.mappings;
    update_indexes();
    return *this;
}

std::vector<std::vector<LongReadMapping>> LongReadsMapper::filter_mappings_by_size_and_id(int64_t size, float id) const {
    std::vector<std::vector<LongReadMapping>> filtered_read_mappings;
    filtered_read_mappings.resize(datastore.size());
    for (auto &m:mappings){
        if (m.nEnd-m.nStart>=size and 100.0*m.score/(m.nEnd-m.nStart)>=id){
            filtered_read_mappings[m.read_id].emplace_back(m);
        }
    }
    return filtered_read_mappings;
}

std::vector<LongReadMapping>
LongReadsMapper::filter_and_chain_matches_by_offset_group(std::vector<LongReadMapping> &matches, bool verbose) const {

    const int32_t OFFSET_BANDWITDH=200;
    if (verbose){
        std::cout<<"Filtering matches:" << std::endl;
        for (const auto &m : matches)
            std::cout << m << std::endl;
    }


    std::vector<match_band> mbo;
    auto node = std::abs(matches[0].node);
    for (const auto &m : matches ){
        mbo.emplace_back(node == m.node, m.qStart-m.nStart, m.qEnd-m.nEnd, m.nEnd-m.nStart, m.score);
    }

    for (auto &m : mbo) {
        auto mino=std::min(m.min_offset,m.max_offset)-OFFSET_BANDWITDH;
        auto maxo=std::max(m.min_offset, m.max_offset)+OFFSET_BANDWITDH;
        m.min_offset = mino;
        m.max_offset = maxo;
    }

    if (verbose){
        std::cout << "Orientations, offsets, lengths and scores:" << std::endl;
        for (const auto &m : mbo){
            std::cout << ((m.dir) ? "true":"false") << ", " << m.min_offset << ", " << m.max_offset << ", " << m.len << ", " << m.score << std::endl;
        }
    }

    auto prev_len = mbo.size()+1;

    while (mbo.size() < prev_len) {
        auto used = std::vector<bool>(mbo.size(), false);
        std::vector<match_band> new_mbo;

        for (int mi=0; mi < mbo.size(); mi++) {
            if (used[mi]){continue;}
            used[mi] = true;
            auto m = mbo[mi];
            for (int mi2 = 0; mi2 < mbo.size(); mi2++) {
                if (used[mi2]){continue;}
                auto m2 = mbo[mi2];
                if (m.min_offset<m2.max_offset and m2.min_offset < m.max_offset){
                    m.min_offset = std::min(m.min_offset,m2.min_offset);
                    m.max_offset = std::max(m.max_offset,m2.max_offset);
                    m.len += m2.len;
                    m.score += m2.score;
                    used[mi2] = true;
                }
            }
            new_mbo.emplace_back(m);
        }
        prev_len = mbo.size();
        mbo = new_mbo;
    }

    std::sort(mbo.begin(), mbo.end(), [](match_band &ma, match_band &mb){return ma.score > mb.score;});

    if (verbose){
        std::cout << "Consolidated orientations, offsets, lengths and scores:\n";
        for (const auto &m : mbo){
            std::cout << ((m.dir) ? "true":"false") << ", " << m.min_offset << ", " << m.max_offset << ", " << m.len << ", " << m.score << std::endl;
        }
    }

    std::vector<match_band> filtered_mbo;
    auto max_score = mbo[0].score;
    std::copy_if(mbo.begin(), mbo.end(), std::back_inserter(filtered_mbo), [max_score](match_band &m){ return m.score >= 0.5f*max_score; });

    if (verbose){
        std::cout << "Filtered orientations, offsets, lengths and scores:\n";
        for (const auto &m : filtered_mbo){
            std::cout << ((m.dir) ? "true":"false") << ", " << m.min_offset << ", " << m.max_offset << ", " << m.len << ", " << m.score << std::endl;
        }
    }

    // Create matches min, max positions for each band that has passed the filter
    std::vector<LongReadMapping> filtered_matches;

    for (const auto &mb : filtered_mbo) {
        LongReadMapping chained;
        chained.nStart=INT32_MAX;
        chained.nEnd=INT32_MIN;
        chained.score=0;
        chained.read_id=matches[0].read_id;
        chained.node=mb.dir ? std::abs(matches[0].node) : -std::abs(matches[0].node);
        for (const auto &m : matches) {
            if ((m.qStart - m.nStart) >= mb.min_offset and (m.qEnd - m.nEnd) <= mb.max_offset) {
                if (chained.nStart>m.nStart){
                    chained.nStart=m.nStart;
                    chained.qStart=m.qStart;
                }
                if (chained.nEnd<m.nEnd){
                    chained.nEnd=m.nEnd;
                    chained.qEnd=m.qEnd;
                }
                chained.score+=m.score;
            }
        }
        filtered_matches.emplace_back(chained);
    }

    if (verbose){
        std::cout<<"Filtered matches:" << std::endl;
        for (const auto &m : filtered_matches)
            std::cout << m << std::endl;
    }

    return filtered_matches;
}

std::vector<LongReadMapping> LongReadsMapper::remove_shadowed_matches(const std::vector<LongReadMapping> &matches) const {
    std::vector<LongReadMapping> filtered_matches;

    for (auto mi1=0;mi1<matches.size();++mi1){
        auto &m1=matches[mi1];
        bool use=true;
        for (auto &m2: matches) {
            if (m2.qStart<=m1.qStart and m2.qEnd>=m1.qEnd and m2.score>m1.score){
                use=false;
                break;
            }
        }
        if (use) filtered_matches.emplace_back(m1);
    }

    return filtered_matches;
}

std::vector<LongReadMapping> LongReadsMapper::improve_read_mappings(const std::vector<LongReadMapping> & input_mappings) const {

    std::vector<LongReadMapping> new_mappings;
    auto imappings = remove_shadowed_matches(input_mappings); // This has to be a copy and can be a bit expensive!
    if (imappings.empty()) return imappings;

    std::map<sgNodeID_t, uint32_t > most_common_map;

    for (const auto &m : imappings) {
        ++most_common_map[std::abs(m.node)];
    }

    std::multimap<uint32_t , sgNodeID_t > most_common = flip_map(most_common_map);

    if (most_common.rbegin()->first == 1) {
        return imappings;
    }

    for (const auto &n : most_common) {
        std::vector<LongReadMapping> nm;
        // This copy will obviously be expensive! It's better to sort the whole structure first and then check if needed to process it.
        std::copy_if(imappings.begin(), imappings.end(), std::back_inserter(nm),
                     [n](const LongReadMapping& m){
                         return n.second == std::abs(m.node);
                     });
        std::sort(nm.begin(), nm.end(),
                  [](const LongReadMapping &a, const LongReadMapping &b) {
                      uint64_t na = std::abs(a.node);
                      uint64_t nb = std::abs(b.node);
                      return std::tie(na, a.qStart) < std::tie(nb, b.qStart);
                  });

        if (nm.size() <= 1){
            new_mappings.insert(new_mappings.end(), nm.begin(), nm.end());
        } else {
            auto fm = filter_and_chain_matches_by_offset_group(nm);
            new_mappings.insert(new_mappings.end(), fm.begin(), fm.end());
        }
    }
    new_mappings=remove_shadowed_matches(new_mappings);
    std::sort(new_mappings.begin(), new_mappings.end(), [](LongReadMapping &a, LongReadMapping &b){
        return a.qStart < b.qStart;
    });

    return new_mappings;
}

std::vector<ReadCacheItem> LongReadsMapper::create_read_paths(const std::vector<sgNodeID_t> &backbone, const std::vector<std::vector<LongReadMapping>> filtered_read_mappings, const ReadPathParams &read_path_params) {
    std::unordered_set<uint64_t> useful_read;
    std::vector<ReadCacheItem> read_cache;
    ReadSequenceBuffer sequenceGetter(datastore);
    for (uint32_t bn = 0; bn < backbone.size(); bn++) {
        for (const auto &read:reads_in_node[std::abs(backbone[bn])]) {
            const auto find_it = useful_read.find(read);
            if (find_it == useful_read.cend()) {
                useful_read.emplace(read);
                std::string seq = sequenceGetter.get_read_sequence(read);
                read_cache.emplace_back(read, seq);
            }

        }
    }
    
    // TODO: HACK! Move this to somewhere it makes more sense
    read_paths.resize(filtered_read_mappings.size());
    
#pragma omp parallel for
    for (uint32_t rcp = 0; rcp < read_cache.size(); rcp++) {
        read_paths[read_cache[rcp].id] = create_read_path(filtered_read_mappings[read_cache[rcp].id], read_path_params, false, read_cache[rcp].seq);
    }

    return read_cache;
}

std::vector<sgNodeID_t> LongReadsMapper::create_read_path(const std::vector<LongReadMapping> mappings, const ReadPathParams &read_path_params, bool verbose, const std::string& read_seq) {
    uint64_t rid=0;
    if (!mappings.empty()) rid=mappings[0].read_id;
    std::vector<sgNodeID_t> read_path;
    if (mappings.empty()) return read_path;

    for (int32_t tid = 0; tid < mappings.size()-1; ++tid) {
        const auto &m1 = mappings[tid];
        const auto &m2 = mappings[tid+1];
        //TODO: Some mappings seem to have not been joined (m1.node == m2.node), needs fixing at mapping stage?
        if (m1.node == m2.node) {continue;}

        read_path.emplace_back(m1.node);

        //this checks the direct connection
        bool need_pathing = true;
        for (const auto &l : sg.get_fw_links(m1.node)) {
            if (l.dest == m2.node) need_pathing = false;
        }

        int32_t ad;
        if (need_pathing){
            ad = (m2.qStart-m2.nStart)-(m1.qEnd+sg.nodes[std::abs(m1.node)].sequence.size()-m1.nEnd);
            auto pd = ad+read_path_params.default_overlap_distance*2; // default_overlap_distance
            auto max_path_size = std::max((int32_t)(pd*read_path_params.path_distance_multiplier),pd+read_path_params.min_path_distance); // path_distance_multiplier, min_path_distance
            if (max_path_size < 0) {
                // TODO: Some mappings have negative distances, needs fixing at mapping stage?

                // XXX: For now simply accept negative distances lower than min_negative_path_distance
                if (std::abs(max_path_size) > read_path_params.min_negative_path_distance) {
                    read_path.emplace_back(0);
                    continue;
                } else {
                    max_path_size = std::abs(max_path_size);
                }
            }
//            max_path_size = 20000;
            if (verbose) std::cout << "\n\nJumping between "<<m1<<" and "<<m2<<std::endl<<"Alignment distance: " << ad << " bp, path distance " << pd << ", looking for paths of up to " << max_path_size << " bp" << std::endl;

            std::unordered_map<std::pair<sgNodeID_t,sgNodeID_t>, std::vector<SequenceGraphPath>>::const_iterator place_in_map;
#pragma omp critical (paths_map)
            {
                place_in_map = all_paths_between.find(std::make_pair(m1.node, m2.node));
            }

            bool found_in_map = place_in_map != all_paths_between.end();
            const std::vector<SequenceGraphPath>& paths = ( found_in_map ?
                    place_in_map->second :
                    sg.find_all_paths_between(m1.node, m2.node, max_path_size, read_path_params.max_path_nodes, false));

            if (!found_in_map) {
                if (verbose) std::cout <<"Adding backbone to collecion!!" <<std::endl;
#pragma omp critical (paths_map)
                {
                    all_paths_between[std::make_pair(m1.node,m2.node)] = paths; // Add path to structure if not there!
                }
            }

            // Max number of paths before giving up
            if (paths.size()>read_path_params.max_distinct_paths) {
                if (verbose) std::cout << "Too many paths" << std::endl;
                read_path.emplace_back(0);
                continue;
            }
            if (verbose)  std::cout << "There are " << paths.size() << " paths" <<std::endl;
            if (paths.empty()){
                read_path.emplace_back(0);
                continue;
            } else {
                //Create a SG with every path as a node
                WorkSpace pws;
                SequenceDistanceGraph psg(pws);
                for (const auto &p : paths) {
                    psg.add_node(Node(p.get_sequence()));
                }

                //Create, parametrise and index a mapper to the paths SG
                auto pm = LongReadsMapper(psg, datastore);
                pm.k = read_path_params.path_mapping_k;
                pm.min_size = std::min(int(ad / 2), read_path_params.max_mapping_path_size); //
                pm.min_chain = read_path_params.min_mapping_chain;
                pm.max_jump = read_path_params.max_mapping_jump;
                pm.max_delta_change = read_path_params.max_mapping_delta_change;
                auto nkindex=NKmerIndex(psg,read_path_params.path_mapping_k,paths.size() * read_path_params.kmer_filter_multiplier + 1);

                //Align the read to the paths SG and pick the best(s) alignement(s)
                std::string seq = read_seq;
                if (read_seq.empty()) {
                    ReadSequenceBuffer sequenceGetter(datastore);
                    seq = sequenceGetter.get_read_sequence(rid);
                }
                auto subs_string = seq.substr(m1.qEnd - 199, m2.qStart + read_path_params.path_mapping_k + 199 - m1.qEnd + 199);
                auto path_mappings = pm.map_sequence(nkindex,subs_string.data());
                std::sort(path_mappings.begin(), path_mappings.end(),
                          [](const LongReadMapping &a, const LongReadMapping &b) { return a.score > b.score; });

                if (path_mappings.empty()) {
                    read_path.emplace_back(0);//paths do not map, add a gap
                    continue;
                } else {
                    if (verbose) {
                        for (int i = 0; i < 10 and i < path_mappings.size(); i++) {
                            if (path_mappings[i].node < 0) continue;
                            std::cout << path_mappings[i] << std::endl;
                            std::copy(paths[path_mappings[i].node - 1].nodes.begin(),
                                      paths[path_mappings[i].node - 1].nodes.end(),
                                      std::ostream_iterator<sgNodeID_t>(std::cout, ", "));
                            std::cout << std::endl;
                        }
                    }


                    if (verbose)std::cout << "Creating winners" << std::endl;
                    std::vector<std::vector<sgNodeID_t>> winners;
                    for (const auto &a: path_mappings) {
                        if (a.score == path_mappings[0].score) {
                            if (a.node < 0) {
                                winners.clear();
                                break;
                            } else {
                                winners.emplace_back(paths[a.node - 1].nodes);
                            }
                        }
                    }

                    if (winners.empty()) {
                        read_path.emplace_back(0);
                        continue;
                    }

                    //Add nodes from winner(s) to path
                    if (winners.size() == 1) {
                        // Single path, all its nodes go to the path
                        if (verbose) std::cout << "Single winner path" << std::endl;
                        read_path.insert(read_path.end(), winners[0].begin(), winners[0].end());
                    } else {
                        if (verbose) std::cout << "Winner paths: " << std::endl;
                        for (const auto &w: winners) {
                            if (verbose) {
                                std::copy(w.cbegin(), w.cend(), std::ostream_iterator<sgNodeID_t>(std::cout, ", "));
                                std::cout << std::endl;
                            }
                        }
                        if (verbose) std::cout << std::endl;
                        std::sort(winners.begin(), winners.end(), [](const std::vector<sgNodeID_t> &a, const std::vector<sgNodeID_t> &b){ return a.size() < b.size(); });
                        std::vector<sgNodeID_t> shared_winner_nodes;
                        for (int pos = 0; pos < winners[0].size(); ++pos) {
                            bool sharing = true;
                            auto &current_node = winners[0][pos];
                            for (const auto &w : winners) {
                                if (w[pos] != current_node) {
                                    sharing=false;
                                    break;
                                }
                            }
                            if (sharing) shared_winner_nodes.emplace_back(current_node);
                            else break;
                        }

                        std::vector<sgNodeID_t> back_shared_winner_nodes;
                        for (int pos = 0; pos < winners[0].size() ; ++pos) {
                            bool sharing = true;
                            auto &current_node = winners[0][winners[0].size()-pos-1];
                            for (const auto &w : winners) {
                                if (w[w.size()-pos-1] != current_node) {
                                    sharing=false;
                                    break;
                                }
                            }
                            if (sharing) back_shared_winner_nodes.emplace_back(current_node);
                            else break;
                        }

                        if (verbose) {
                            std::cout << std::endl << "Shared nodes in winner paths: " << std::endl;
                        }
                        read_path.insert(read_path.end(), shared_winner_nodes.begin(), shared_winner_nodes.end());
                        if (verbose) {
                            std::copy(shared_winner_nodes.cbegin(), shared_winner_nodes.cend(), std::ostream_iterator<sgNodeID_t>(std::cout, ", "));
                            std::cout << "0, ";
                        }
                        read_path.emplace_back(0);//as of now, if multiple winners, just put a full gap
                        read_path.insert(read_path.end(), back_shared_winner_nodes.rbegin(), back_shared_winner_nodes.rend());
                        if (verbose) {
                            std::copy(back_shared_winner_nodes.crbegin(), back_shared_winner_nodes.crend(), std::ostream_iterator<sgNodeID_t>(std::cout, ", "));
                            std::cout << std::endl;
                        }
                    }
                }
            }
        } else {
            if (verbose) std::cout << "Nodes directly connected, no pathing required\n";
        }
    }
    read_path.emplace_back(mappings.back().node);
    return read_path;
}

LongReadsMapper::LongReadsMapper(const LongReadsDatastore &ds, const LongReadsMapper &o) :
datastore(ds),
sg(o.sg),
read_paths(o.read_paths),
first_mapping(o.first_mapping),
k(o.k)
{}

std::ostream &operator<<(std::ostream &os, const LongReadsMapper &lorm) {
    os << "LongReadsMapper" << std::endl;
    lorm.print_status();
}
