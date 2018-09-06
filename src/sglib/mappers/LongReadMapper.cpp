//
// Created by Luis Yanes (EI) on 22/03/2018.
//

#include <cmath>
#include <sglib/mappers/LongReadMapper.hpp>
#include <sglib/logger/OutputLog.h>

#include <sglib/utilities/omp_safe.hpp>
#include <atomic>

const bsgVersion_t LongReadMapper::min_compat = 0x0001;

void LongReadMapper::map_reads(std::unordered_set<uint32_t> readIDs) {
    if (assembly_kmers.empty()) assembly_kmers.generate_index(sg);
    std::vector<std::vector<LongReadMapping>> thread_mappings(omp_get_max_threads());
    std::atomic<uint32_t > num_reads_done(0);
#pragma omp parallel
    {
        StringKMerFactory skf(k);
        std::vector<std::pair<bool, uint64_t>> read_kmers;

        auto & private_results=thread_mappings[omp_get_thread_num()];
        BufferedSequenceGetter sequenceGetter(datastore);
        std::vector<std::vector<std::pair<uint32_t, bool>>> node_matches;
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
            bool print_window(false);

            const auto query_sequence(sequenceGetter.get_read_sequence(readID));
            if ( query_sequence.size()< 4*window_size) {
                continue;
            }

            std::map<uint32_t , uint32_t > node_score;
            //===== Create a vector of nodes for every k-mer in the read.
            read_kmers.clear();
            skf.create_kmers(query_sequence, read_kmers); //XXX: Ns will produce "deletion-windows"

            if (node_matches.size() < read_kmers.size()) {
                node_matches.resize(read_kmers.size());
            }

            for (auto i=0;i<read_kmers.size();++i){
                node_matches[i].clear();
                auto first = assembly_kmers.find(read_kmers[i].second);
                for (auto it = first; it != assembly_kmers.end() && it->kmer == read_kmers[i].second; ++it) {
                    if (read_kmers[i].first != std::signbit(it->offset)){
                        node_matches[i].emplace_back(it->contigID, true);
                    }
                    else {
                        node_matches[i].emplace_back(it->contigID, false);
                    }
                    ++node_score[it->contigID];
                }
            }



            //===== Rank nodes by how many kmer-hits they have on the read, choose the top N
            std::vector<std::pair<uint32_t , uint32_t>> score_node;
            for (auto &ns:node_score){
                score_node.emplace_back(ns.second, ns.first);
            }
            std::sort(score_node.rbegin(),score_node.rend());

            if (score_node.size()>max_num_score_nodes) {
                score_node.resize(max_num_score_nodes);
            }

            //===== For every window (sliding by win_size/2) save total of kmers in every node of the top N
//            bool print_csv(false);
//            std::ofstream readout("read" + std::to_string(readID+1) + ".csv");
//        if (print_csv) {
//            readout << "window#";
//            for (auto n:score_node) readout << "," << "seq" <<n.second;
//            readout << std::endl;
//        }

            std::vector< std::pair<int32_t, float> > winners;
            for (int wp_i = 0 , w_i = 0; wp_i < read_kmers.size(); ++w_i) {
                std::map<uint32_t,float> win_node_score;
                auto win_start(w_i * window_slide);
                auto win_end(w_i * window_slide + window_size);

                for (wp_i = win_start; wp_i < win_end and wp_i < read_kmers.size(); wp_i++) {
                    for (auto matchnode:node_matches[wp_i]) {
                        if (matchnode.second) win_node_score[matchnode.first] += 1.f/node_matches[wp_i].size();
                        else win_node_score[-matchnode.first] += 1.f/node_matches[wp_i].size();
                    }
                }

                std::vector<std::pair<int32_t,float>> win_node_vec(win_node_score.begin(),win_node_score.end());

                std::sort(win_node_vec.begin(), win_node_vec.end(),
                          [](const std::pair<uint32_t ,float> &a, const std::pair<uint32_t ,float> &b) {return std::abs(a.second) > std::abs(b.second);});

                std::pair<int32_t, float> winner(0,0.0f);

                if (print_window) sglib::OutputLog() << "Window " << w_i+1 << " " << win_start << " to " << win_end << " ";
//            if (print_csv) {
//                readout << w_i;
//            }
                int lt_minkmers(0);
                int spread_bellow_min(0);
                if (win_node_vec.size()>0 and std::abs(win_node_vec[0].second) > min_kmers) {
                    if (win_node_vec.size()<2 or std::abs(win_node_vec[0].second*min_spread) > std::abs(win_node_vec[1].second)) {
                        if (print_window) std::cout << win_node_vec[0].first-1 << " " << win_node_vec[0].second;
                        winner = win_node_vec[0];
                    } else {
                        spread_bellow_min++;
                    }
                } else {
                    lt_minkmers++;
                }
                if (print_window) sglib::OutputLog() << " lt_minkmers " << lt_minkmers << " spread " << spread_bellow_min << std::endl;

                winners.emplace_back(winner);
//                if (print_csv) {
//                    for (auto n:score_node) {
//                        readout << "," << win_node_score[n.second];
//                    }
//                    readout << std::endl;
//                }
            }

            if (winners.empty()) {
                continue;
            }
            int rp(0),wp(0),last_winner(0);
            while(rp < winners.size() and winners[rp].first==0){rp++;} // Find first non zero
            if (rp == winners.size()) continue;
            wp=rp;
            last_winner=winners[wp].first;
            wp++;
            rp++;
            for (; rp < winners.size(); rp++) {
                while (rp < winners.size() and winners[wp].first == winners[rp].first) {
                    rp++;
                }
                if (rp-wp > 1) {
                    if (winners[wp].first != 0 and winners[wp].first != last_winner) {
                        last_winner = winners[wp].first;
                        private_results.emplace_back(winners[wp].first, readID, 0, 0, 0, 0); //TODO: Use the correct positions!
                    }
                }
                wp=rp;
            }
        }
    }
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

void LongReadMapper::update_graph_index() { assembly_kmers.generate_index(sg); }
