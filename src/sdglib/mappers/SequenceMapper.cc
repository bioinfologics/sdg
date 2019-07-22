//
// Created by Bernardo Clavijo (EI) on 2019-07-12.
//

#include "SequenceMapper.hpp"

SequenceMapper::SequenceMapper(const DistanceGraph &_dg,uint8_t _k, int _min_size, int _min_chain, int _max_jump, int _max_delta_change):sg(_dg) {
    k=_k;
    min_size=_min_size;
    min_chain=_min_chain;
    max_jump=_max_jump;
    max_delta_change=_max_delta_change;
    update_graph_index();
}

std::vector<LongReadMapping> SequenceMapper::map_sequence(const char * query_sequence_ptr, sgNodeID_t seq_id) {
    std::vector<LongReadMapping> private_results;
    StreamKmerFactory skf(k);
    std::vector<std::pair<bool, uint64_t>> read_kmers;
    std::vector<std::vector<std::pair<int32_t, int32_t>>> node_matches; //node, offset
    read_kmers.clear();
    std::vector<unsigned char> candidate_counts(sg.sdg.nodes.size()*2);
    skf.produce_all_kmers(query_sequence_ptr, read_kmers);
    if (sat_kmer_index) {
        get_sat_kmer_matches(node_matches, read_kmers);
    } else {
        get_all_kmer_matches(node_matches,read_kmers);
    }

    //========== 2. Find match candidates in fixed windows ==========
    count_candidates(candidate_counts, node_matches,read_kmers.size());

    //========== 3. Create alignment blocks from candidates ==========
    return alignment_blocks(seq_id,node_matches,read_kmers.size(), candidate_counts);

}

void SequenceMapper::update_graph_index(int filter_limit, bool verbose) {
    if (sat_kmer_index) {
        if (verbose) std::cout<<"updating satindex with k="<<std::to_string(k)<<std::endl;
        sat_assembly_kmers=SatKmerIndex(k);
        sat_assembly_kmers.generate_index(sg.sdg, filter_limit, verbose);
    } else {
        if (verbose) std::cout<<"updating nkindex with k="<<std::to_string(k)<<std::endl;
        assembly_kmers = NKmerIndex(k);
        assembly_kmers.generate_index(sg.sdg, filter_limit, verbose);
    }
}

/*std::vector<LongReadMapping>
SequenceMapper::filter_and_chain_matches_by_offset_group(std::vector<LongReadMapping> &matches, bool verbose) {

    const int32_t OFFSET_BANDWITDH=200;
    if (verbose){
        std::cout<<"Filtering matches:" << std::endl;
        for (const auto &m : matches)
            std::cout << m << std::endl;
    }


    std::vector<struct match_band> mbo;
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
}*/

std::vector<LongReadMapping> SequenceMapper::remove_shadowed_matches(std::vector<LongReadMapping> &matches,
                                                                      bool verbose) {
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

/*std::vector<LongReadMapping> SequenceMapper::improve_read_filtered_mappings(uint32_t rid, bool correct_on_ws) {

    std::vector<LongReadMapping> new_mappings;
    auto mappings = remove_shadowed_matches(filtered_read_mappings[rid]); // This has to be a copy and can be a bit expensive!
    if (mappings.empty()) return mappings;

    std::map<sgNodeID_t, uint32_t > most_common_map;

    for (const auto &m : mappings) {
        ++most_common_map[std::abs(m.node)];
    }

    std::multimap<uint32_t , sgNodeID_t > most_common = flip_map(most_common_map);

    if (most_common.rbegin()->first == 1) {
        filtered_read_mappings[rid] = mappings;
        return mappings;
    }

    for (const auto &n : most_common) {
        std::vector<LongReadMapping> nm;
        // This copy will obviously be expensive! It's better to sort the whole structure first and then check if needed to process it.
        std::copy_if(mappings.begin(), mappings.end(), std::back_inserter(nm),
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

    if (correct_on_ws) {
        filtered_read_mappings[rid] = new_mappings;
    }
    return new_mappings;
}

void SequenceMapper::filter_mappings_by_size_and_id(int64_t size, float id) {
    filtered_read_mappings.clear();
    filtered_read_mappings.resize(datastore.size());
    for (auto &m:mappings){
        if (m.nEnd-m.nStart>=size and 100.0*m.score/(m.nEnd-m.nStart)>=id){
            filtered_read_mappings[m.read_id].emplace_back(m);
        }
    }
}*/

void SequenceMapper::get_sat_kmer_matches(std::vector<std::vector<std::pair<int32_t, int32_t>>> &matches, std::vector<std::pair<bool, uint64_t>> &read_kmers) {
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
                offset=( (int) sg.sdg.nodes[std::llabs(sat_assembly_kmers.contig_offsets[it].contigID)].sequence.size() ) - std::abs(offset);
            }
            else offset=std::abs(offset)-1;
            matches[i].emplace_back(node, offset);
        }
        if (matches[i].empty()) ++no_match; //DEBUG
        else if (matches[i].size()==1) ++single_match; //DEBUG
        else ++multi_match; //DEBUG
    }
}

void SequenceMapper::get_all_kmer_matches(std::vector<std::vector<std::pair<int32_t, int32_t>>> & matches, std::vector<std::pair<bool, uint64_t>> & read_kmers) {
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
                offset=( (int) sg.sdg.nodes[std::llabs(it->contigID)].sequence.size() ) - std::abs(offset);
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

void SequenceMapper::count_candidates(std::vector<unsigned char> &candidate_counts,
                                       std::vector<std::vector<std::pair<int32_t, int32_t>>> &matches,
                                       uint32_t read_kmers_size){
    //beware you need the size! otherwise there is a size from the matches and not from the reads!
    for (auto i=0;i<read_kmers_size;++i) {
        for (auto m:matches[i]) {
            uint32_t linear_node = (m.first>0) ? m.first : sg.sdg.nodes.size()+std::abs(m.first);
            uint16_t count_so_far=candidate_counts[ linear_node ]+1;
            candidate_counts[ linear_node ] = (count_so_far > UINT8_MAX) ? UINT8_MAX : count_so_far;
        }
    }
}

std::vector<LongReadMapping> SequenceMapper::alignment_blocks(uint32_t readID,
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
            uint32_t linear_node = (m.first>0) ? m.first : sg.sdg.nodes.size()+std::abs(m.first);
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