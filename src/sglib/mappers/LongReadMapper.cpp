//
// Created by Luis Yanes (EI) on 22/03/2018.
//

#include <cmath>
#include <sglib/mappers/LongReadMapper.hpp>
#include <sglib/logger/OutputLog.h>

#include <sglib/utilities/omp_safe.hpp>
#include <atomic>

const bsgVersion_t LongReadMapper::min_compat = 0x0001;

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
        StringKMerFactory skf(k);
        std::vector<std::pair<bool, uint64_t>> read_kmers;

        auto & private_results=thread_mappings[omp_get_thread_num()];
        BufferedSequenceGetter sequenceGetter(datastore);
        std::vector<std::vector<std::pair<int32_t, int32_t>>> node_matches; //node, offset
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
                (readIDs.empty() and !read_to_mappings[readID].empty()/*No selection and read is mapped*/)) {
                continue;
            }

            const auto query_sequence(sequenceGetter.get_read_sequence(readID));
            if ( query_sequence.size()< 4*window_size) {
                continue;
            }
            //std::cout<<"processing read #"<<readID<<std::endl;

            //get the kmers from the read
            //std::cout<<"creating kmers"<<std::endl;
            read_kmers.clear();
            skf.create_kmers(query_sequence, read_kmers); //XXX: Ns will produce "deletion-windows" in the read
            //std::cout<<"kmers created"<<std::endl;
            //populate vector of vectors with the matches for every kmer
            if (node_matches.size() < read_kmers.size()) {
                node_matches.resize(read_kmers.size());
            }

            //std::cout<<"read has "<<read_kmers.size()<<" kmers"<<std::endl;
            uint32_t total_windows=(read_kmers.size()-window_size)/window_slide+1;
            //std::cout<<"read has "<<total_windows<<" windows of "<<window_size<<" bp, every "<<window_slide<<" bp"<<std::endl;
            std::map<uint32_t , uint32_t > node_score[total_windows];
            //std::cout<<"adding hits and scores"<<std::endl;
            for (auto i=0;i<read_kmers.size();++i){
                node_matches[i].clear();
                auto first = assembly_kmers.find(read_kmers[i].second);
                //std::cout<<"Processing hits for position #"<<i<<" read k-mer orientation "<<(read_kmers[i].first ? "FW":"RC")<<std::endl;
                for (auto it = first; it != assembly_kmers.end() && it->kmer == read_kmers[i].second; ++it) {
                    //XXX: yes, offset 0 is screwed up because -0 can't be used for reverse, so it is not used
                    //std::cout<<" hit to "<<it->contigID<<" : "<<it->offset<<std::endl;
                    int32_t offset=it->offset; //so far, this is +1 and the sign indicate direction of kmer in contig
                    sgNodeID_t node=it->contigID; //so far, this is always positive
                    if (read_kmers[i].first != (offset>0) ) {
                        node=-node;
                        offset=( (int) sg.nodes[std::llabs(it->contigID)].sequence.size() ) - std::abs(offset);
                    }
                    else offset=std::abs(offset)-1;
                    //std::cout<<" node "<<node<<", offset "<<offset<<std::endl;
                    node_matches[i].emplace_back(node, offset);
                    //std::cout<<"  Adding result to windows "<<std::max(0,(i-window_size)/window_slide)<< " to "<<std::min(i/window_slide,(int)total_windows-1)<<std::endl;
                    for (auto win=std::max(0,(i-window_size)/window_slide+1);win<std::min(i/window_slide,(int)total_windows-1);++win) {
                        win=++node_score[win][node];
                    }
                }
            }


            //std::cout<<"picking candidates"<<std::endl;
            //===== Rank nodes by how many kmer-hits they have on the read, choose the top N
            //pick winners in all windows
            std::map<uint32_t , int32_t > node_total_score;
            for (auto win=0;win<total_windows;++win){
                std::vector<std::pair<uint32_t , int32_t>> winscore_node;
                for (auto &ns:node_score[win]){
                    winscore_node.emplace_back(ns.second, ns.first);
                }
                std::sort(winscore_node.rbegin(),winscore_node.rend());

                for (auto i=0;i<std::min((int)winscore_node.size(),10);++i)
                    ++node_total_score[winscore_node[i].second];
            }
            //std::vector<std::pair<uint32_t , int32_t>> score_node;
            std::unordered_set<int32_t> candidates;
            for (auto &ns:node_total_score){
                if (ns.second>3) candidates.insert(ns.first);
            }

            /*if (score_node.size()>max_num_score_nodes) {
                score_node.resize(max_num_score_nodes);
            }*/

            //===== For every window (sliding by win_size/2) score
            //std::cout<<"picking window winners among "<<candidates.size()<<" candidates:";
            //for (auto c:candidates) std::cout<<" "<<c;
            //std::cout<<std::endl;

            std::vector< std::pair<int32_t, float> > winners;

            for (int window = 0; window<total_windows; ++window) {
                std::map<int32_t,float> win_node_score;

                auto win_start(window * window_slide);
                auto win_end(std::min((int)read_kmers.size(),window * window_slide + window_size));

                //aggregate node scores of 1/kmer-hits for every position, win_node_score keys are directional nodes
                for (auto wp_i = win_start; wp_i < win_end ; ++wp_i) {
                    int total_candidate_matches=0;
                    for (auto matchnode:node_matches[wp_i])
                        if (candidates.count(matchnode.second > 0 ? matchnode.first:-matchnode.first)) ++total_candidate_matches;
                    for (auto matchnode:node_matches[wp_i]) {
                        if (candidates.count(matchnode.second>0 ? matchnode.first:-matchnode.first))
                            win_node_score[matchnode.first] += 1.f/(std::pow(2,total_candidate_matches));
                    }
                }
                //TODO: no need to do a vector, sort, etc only to consider second-best score!

                std::vector<std::pair<int32_t,float>> win_node_vec(win_node_score.begin(),win_node_score.end());

                std::sort(win_node_vec.begin(), win_node_vec.end(),
                          [](const std::pair<uint32_t ,float> &a, const std::pair<uint32_t ,float> &b) {return std::abs(a.second) > std::abs(b.second);});

                std::pair<int32_t, float> winner(0,0.0f);

                if (win_node_vec.size()>0) {
                    if (win_node_vec[0].second < min_score) {
                        ++window_low_score;
                        winner.second=1;
                    }
                    else if (win_node_vec.size() > 1 and
                             win_node_vec[1].second * 100 >= win_node_vec[0].second * second_best_score_pct) {
                        ++window_close_second;
                        winner.second=2;
                    }
                    else {
                        winner = win_node_vec[0];
                        ++window_hit;
                    }
                }
                //std::cout<<"winner: "<<winner.first<<std::endl;
                winners.emplace_back(winner);
            }
            if (!detailed_log.empty()){
#pragma omp critical
                {
                    dl << readID << ": ";
                    for (auto &w:winners) dl<<" "<<w.first<<"("<<w.second<<")";
                    dl << std::endl;
                }
            }
            //std::cout<<"winners picked, writing mappings"<<std::endl;
            //Now transform windows of winners into mappings
            {
                //std::cout<<"Aggregating among "<<winners.size()<<" winners:";
                //for (auto &w:winners)std::cout<<"  "<<w.first<< " ("<<w.second<<")";
                //std::cout<<std::endl;
                auto prs=private_results.size();
                auto win = 0;
                for (; win < winners.size(); ++win) {
                    auto win_start = win;
                    while (win_start < winners.size() and winners[win_start].first == 0) ++win_start;
                    if (win_start == winners.size()) break;
                    //std::cout<<" extending window from: " << win_start <<": "<<winners[win_start].first<<std::endl;
                    win=win_start;
                    auto win_end=win;
                    int matches=0;
                    while (win < winners.size() and
                           (winners[win].first == 0 or winners[win].first == winners[win_start].first)) {
                        if (winners[win].first!=0) {
                            win_end=win;
                            ++matches;
                        }
                        ++win;
                    }
                    //std::cout<<"  window extended to: " << win_end <<": "<<winners[win_end].first<<std::endl;
                    //Ok, so now every winner in winners[win_start:win_end] is either 0 or the winner node;
                    //TODO: use median or >20% ^ <80% to choose start/finish
                    if (matches>5) private_results.emplace_back(winners[win_start].first, readID, 0, 0,
                                                 window_slide*win_start, window_slide*win_end+window_size);//TODO: Use the correct positions: first kmer from first, last kmer fomr last
                }
                if (prs==private_results.size()) ++no_matches;
                else if (prs+1==private_results.size()) ++single_matches;
                else ++multi_matches;
            }
        }
    }
    //TODO: report read and win coverage by winners vs. by loosers
    sglib::OutputLog()<<"Read window results:    "<<window_low_score<<" low score    "<<window_close_second<<" close second    "<<window_hit<<" hits"<<std::endl;
    sglib::OutputLog()<<"Read results:    "<<no_matches<<" no match    "<<single_matches<<" single match    "<<multi_matches<<" multi matches"<<std::endl;
    for (auto & threadm : thread_mappings) {
        mappings.insert(mappings.end(),threadm.begin(),threadm.end());
    }
    sglib::OutputLog() << "Updating mapping indexes" <<std::endl;
    update_indexes_from_mappings();
}

void LongReadMapper::update_indexes_from_mappings() {
    reads_in_node.clear();
    reads_in_node.resize(sg.nodes.size());
    read_to_mappings.clear();
    read_to_mappings.resize(datastore.size());
    for (auto i=0;i<mappings.size();++i) {
        auto &m=mappings[i];
        read_to_mappings[m.read_id].push_back(i);
        reads_in_node[std::abs(m.node)].push_back(m.read_id);
    }

    sglib::OutputLog() << "Finished with " << mappings.size() << " mappings" << std::endl;
}

LongReadMapper::LongReadMapper(SequenceGraph &sg, LongReadsDatastore &ds, uint8_t k)
        : sg(sg), k(k), datastore(ds), assembly_kmers(k) {

    reads_in_node.resize(sg.nodes.size());
    read_to_mappings.resize(datastore.size());
}

LongReadMapper::~LongReadMapper() {}

void LongReadMapper::read(std::string filename) {
    // Read the mappings from file
    sglib::OutputLog() << "Reading long read mappings" << std::endl;
    std::ifstream inf(filename, std::ios_base::binary);
    auto mapSize(mappings.size());
    inf.read(reinterpret_cast<char *>(&k), sizeof(k));
    inf.read(reinterpret_cast<char *>(&mapSize), sizeof(mapSize));
    mappings.reserve(mapSize);
    inf.read(reinterpret_cast<char*>(mappings.data()), mappings.size()*sizeof(LongReadMapping));

    sglib::OutputLog() << "Updating read mapping indexes!" << std::endl;
    update_indexes_from_mappings();
    sglib::OutputLog() << "Done!" << std::endl;
}

void LongReadMapper::read(std::ifstream &inf) {
    auto mapSize(mappings.size());
    inf.read(reinterpret_cast<char *>(&k), sizeof(k));
    inf.read(reinterpret_cast<char *>(&mapSize), sizeof(mapSize));
    mappings.resize(mapSize);
    inf.read(reinterpret_cast<char*>(mappings.data()), mappings.size()*sizeof(LongReadMapping));

    sglib::OutputLog() << "Updating read mapping indexes!" << std::endl;
    update_indexes_from_mappings();
    sglib::OutputLog() << "Done!" << std::endl;
}

void LongReadMapper::write(std::string filename) {
    // Write mappings to file
    sglib::OutputLog() << "Dumping long read mappings" << std::endl;
    std::ofstream outf(filename, std::ios_base::binary);
    auto mapSize(mappings.size());
    outf.write(reinterpret_cast<const char *>(&k), sizeof(k));
    outf.write(reinterpret_cast<const char *>(&mapSize), sizeof(mapSize));
    outf.write(reinterpret_cast<const char*>(mappings.data()), mappings.size()*sizeof(LongReadMapping));
    sglib::OutputLog() << "Done!" << std::endl;
}

void LongReadMapper::write(std::ofstream &ofs) {
    auto mapSize(mappings.size());
    ofs.write(reinterpret_cast<char *>(&k), sizeof(k));
    ofs.write(reinterpret_cast<char *>(&mapSize), sizeof(mapSize));
    mappings.reserve(mapSize);
    ofs.write(reinterpret_cast<char*>(mappings.data()), mappings.size()*sizeof(LongReadMapping));
    sglib::OutputLog() << "Done!" << std::endl;
}

LongReadMapper LongReadMapper::operator=(const LongReadMapper &other) {
    if (&sg != &other.sg and &datastore != &other.datastore) { throw ("Can only copy paths from the same SequenceGraph"); }
    if (&other == this) {
        return *this;
    }
    mappings = other.mappings;
    update_indexes_from_mappings();
    return *this;
}

std::vector<uint64_t> LongReadMapper::get_node_read_ids(sgNodeID_t nodeID) const {
    std::vector<uint64_t> rpin;
    nodeID=llabs(nodeID);
    rpin.reserve(reads_in_node[nodeID].size());
    for (auto &rm:reads_in_node[nodeID]){
        rpin.emplace_back(rm);
    }
    std::sort(rpin.begin(),rpin.end());
    return rpin;
}

void LongReadMapper::update_graph_index() {
    std::cout<<"updating index with k="<<std::to_string(k)<<std::endl;
    assembly_kmers=NKmerIndex(k);
    assembly_kmers.generate_index(sg);
}
