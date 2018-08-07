//
// Created by Luis Yanes (EI) on 22/03/2018.
//

#include <cmath>
#include <sglib/mappers/LongReadMapper.hpp>

#ifdef _OPENMP
#include <omp.h>
#include <parallel/algorithm>
#include <sglib/logger/OutputLog.h>

#endif

void LongReadMapper::map_reads(std::unordered_set<uint32_t> readIDs) {
    std::vector<std::vector<LongReadMapping>> thread_mappings(omp_get_max_threads());

#pragma omp parallel
    {
        StringKMerFactory skf(k);
        std::vector<std::pair<bool, uint64_t>> read_kmers;
        std::ofstream mapping_output( "mappings" + std::to_string(omp_get_thread_num()) + ".tsv");
        mapping_output << "query\tnodes\tscores" << std::endl;

#pragma omp for
        for (uint32_t readID = 1; readID < datastore.size(); ++readID) {
            if (!((readIDs.size()>0 and readIDs.count(readID)>0) or (readIDs.empty() and read_to_mappings[readID].empty())))
                continue;

            bool print_csv(false);
            bool print_window(false);

            const auto &line(datastore.get_read_sequence(readID));
            if ( line.size()< 4*window_size) continue;
            sglib::OutputLog() << "Processing read " << readID << " " << line.size() << std::endl;
            std::map<uint32_t , uint32_t > node_score;
            //===== Create a vector of nodes for every k-mer in the read.
            read_kmers.clear();
            skf.create_kmers(line, read_kmers); //XXX: Ns will produce "deletion-windows"
            sglib::OutputLog() <<"Generating node matches... SLOOOOOOW"<<std::endl;
            std::vector<std::vector<std::pair<uint32_t, bool>>> node_matches(read_kmers.size());
            for (auto i=0;i<read_kmers.size();++i){
                auto first = std::lower_bound(assembly_kmers.begin(), assembly_kmers.end(), read_kmers[i].second, KmerIDX::ltKmer());
                for (auto it = first; it != assembly_kmers.end() && it->kmer == read_kmers[i].second; ++it) {
                    if (read_kmers[i].first != std::signbit(it->pos)) node_matches[i].emplace_back(it->contigID, true);
                    else node_matches[i].emplace_back(it->contigID, false);
                    ++node_score[it->contigID];
                }
            }
            sglib::OutputLog()<<"DONE"<<std::endl;



            //===== Rank nodes by how many kmer-hits they have on the read, choose the top N
            sglib::OutputLog()<<"Scoring nodes"<<std::endl;
            std::vector<std::pair<uint32_t , uint32_t>> score_node;
            for (auto &ns:node_score){
                score_node.emplace_back(ns.second, ns.first);
            }
            std::sort(score_node.rbegin(),score_node.rend());

            if (score_node.size()>max_num_score_nodes){
                sglib::OutputLog() << "This read contains kmers from " << score_node.size() << " nodes, limiting to " << max_num_score_nodes << std::endl;
                score_node.resize(max_num_score_nodes);
            } else {
                sglib::OutputLog() << "This read contains kmers from " << score_node.size() << " nodes" << std::endl;
            }
            sglib::OutputLog()<<"DONE"<<std::endl;


            //===== For every window (sliding by win_size/2) save total of kmers in every node of the top N
            sglib::OutputLog()<<"Collecting votes"<<std::endl;
//        if (print_csv) {
//            std::ofstream readout("read" + std::to_string(i / 2) + ".csv");
//            readout << "window#";
//            for (auto n:score_node) readout << "," << node_names[(n.second) / 2];
//            readout << std::endl;
//        }

            mapping_output << readID << "\t";
            std::vector< std::pair<uint32_t, float> > winners;
            for (int wp_i = 0 , w_i = 0; wp_i < read_kmers.size(); ++w_i) {
                std::map<uint32_t,float> win_node_score;
                auto win_start(w_i * window_slide);
                auto win_end(w_i * window_slide + window_size);

                for (wp_i = win_start; wp_i < win_end and wp_i < read_kmers.size(); wp_i++) {
                    for (auto matchnode:node_matches[wp_i]) {
                        if (matchnode.second) win_node_score[matchnode.first] += 1.f/node_matches[wp_i].size();
                        else win_node_score[matchnode.first] -= 1.f/node_matches[wp_i].size();
                    }
                }

                std::vector<std::pair<int,float>> win_node_vec(win_node_score.begin(),win_node_score.end());

                std::sort(win_node_vec.begin(), win_node_vec.end(),
                          [](const std::pair<uint32_t ,float> &a, const std::pair<uint32_t ,float> &b) {return std::abs(a.second) > std::abs(b.second);});

                std::pair<uint32_t, float> winner(0,0.0f);

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
//            if (print_csv) {
//                for (auto n:score_node) {
//                    readout << "," << win_node_score[n.second];
//                }
//                readout << std::endl;
//            }
            }
            sglib::OutputLog()<<"DONE"<<std::endl;
            mapping_output<<( (winners[0].second>=0) ? "" : "-" )<< ( (winners[0].first == 0 ) ? "0": sg.oldnames[(winners[0].first)/2]);
            for (int i = 1; i < winners.size(); ++i) {
                mapping_output<<","<<( (winners[i].second>=0) ? "" : "-" )<< ( (winners[i].first == 0 ) ? "0": sg.oldnames[(winners[i].first)/2]);
            }
            mapping_output << "\t";
            mapping_output<<winners[0].second;
            for (int i = 1; i < winners.size(); ++i) {
                mapping_output<< ","<<winners[i].second;
            }

            int rp(0),wp(0),last_winner(0);
            sglib::OutputLog() << "Path " << readID << " ";
            while(winners[rp].first==0){rp++;} // Find first non zero
            wp=rp;
            if (rp < winners.size()) std::cout << sg.oldnames[winners[wp].first/2];
            last_winner=winners[wp].first;
            wp++;
            rp++;
            for (; rp < winners.size(); rp++) {
                while (winners[wp].first == winners[rp].first and rp < winners.size()) {
                    rp++;
                }
                if (rp-wp > 1) {
                    if (winners[wp].first != 0 and winners[wp].first != last_winner) {
                        std::cout << "," << ((winners[wp].second>0)?"":"-") << sg.oldnames[winners[wp].first/2];
                        last_winner = winners[wp].first;
                    }
                }
                wp=rp;
            }
            std::cout << std::endl;

            mapping_output << std::endl;
        }
    }
    for (int thread = 0; thread<omp_get_max_threads(); thread++) {
        mappings.reserve(thread_mappings.size());
        for (const auto &mapping : thread_mappings[thread]) {
            mappings.emplace_back(mapping);
        }
    }

    update_indexes_from_mappings();
}

void LongReadMapper::update_indexes_from_mappings() {
    for (std::vector<LongReadMapping>::const_iterator mappingItr = mappings.cbegin(); mappingItr != mappings.cend(); ++mappingItr) {
        auto index = (unsigned long)std::distance(mappings.cbegin(), mappingItr);
        read_to_mappings[mappingItr->read_id].push_back(index);
        mappings_in_node[std::abs(mappingItr->node)].push_back(index);
    }
}

LongReadMapper::LongReadMapper(SequenceGraph &sg, LongReadsDatastore &ds, uint8_t k)
        : sg(sg), k(k), datastore(ds) {

    update_graph_index();

    mappings_in_node.resize(sg.nodes.size());
    read_to_mappings.resize(datastore.size());
}

void LongReadMapper::update_graph_index() {
    assembly_kmers.reserve(110000000);

    StringKMerFactory skf(k);
    std::vector<std::pair<bool,uint64_t > > contig_kmers;

    for (sgNodeID_t n = 1; n < sg.nodes.size(); ++n) {
        if (sg.nodes[n].sequence.size() >= k) {
            contig_kmers.clear();
            skf.create_kmers(sg.nodes[n].sequence, contig_kmers);
            int k_i(0);
            for (const auto &kmer:contig_kmers) {
                assembly_kmers.emplace_back(kmer.second, n, k_i, 1);
                k_i++;
            }
        }
    }
#ifdef _OPENMP
    __gnu_parallel::sort(assembly_kmers.begin(),assembly_kmers.end(), KmerIDX::byKmerContigOffset());
#else
    std::sort(assembly_kmers.begin(),assembly_kmers.end(), kmerPos::byKmerContigOffset());
#endif
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
