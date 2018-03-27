//
// Created by Luis Yanes (EI) on 22/03/2018.
//

#include <sglib/mappers/LongReadMapper.hpp>
#include <sglib/utilities/omp_safe.hpp>
/**
 * Receives a filename containing long reads in FASTQ format, for each read generates a _valid_ (the first path found)
 * path between matched nodes if no path can be found starts a new path, each read can generate multiple paths of
 * *one* or more nodes.
 * @param filename FASTQ file containing long reads
 * @return A vector of SequenceGraphPath containing paths between connected nodes
 */
void LongReadMapper::map_reads(std::unordered_set<uint32_t> readIDs) {

    std::vector< std::vector<LongReadMapping> > read_mappings(datastore.size());
    if (not readIDs.empty())
        sglib::OutputLog()<<readIDs.size()<<" selected reads / "<<datastore.size()-1<<" total"<<std::endl;

#pragma omp parallel
    {
        std::ofstream matchOutput(std::string("thread_")+std::to_string(omp_get_thread_num())+std::string(".match"));
        std::ofstream blockOutput(std::string("thread_")+std::to_string(omp_get_thread_num())+std::string(".block"));
#pragma omp for
        for (uint32_t readID=1;readID<datastore.size();++readID) {
            if ( (!readIDs.empty() and readIDs.count(readID)) or (readIDs.empty() and read_to_mappings[readID].empty()) ) {
                const auto currentReadMappings(map_read
                                                       (readID,
                                                        datastore.get_read_sequence(readID),
                                                        matchOutput,
                                                        blockOutput)
                );
                read_mappings[readID] = currentReadMappings;
            }
        }
    }
    // Initialise reads_in_node, read_to_nodes from results
    std::vector<LongReadMapping>::size_type totalMappings(0);
    for (const auto &rm : read_mappings) {
        totalMappings += rm.size();
    }
    mappings.reserve(totalMappings);
    for (const auto &rm : read_mappings) {
        mappings.insert(mappings.end(), rm.cbegin(), rm.cend());
    }

    for (std::vector<LongReadMapping>::size_type i=1; i < mappings.size(); ++i) {
        mappings_in_node[std::abs(mappings[i].node)].push_back(i);
        read_to_mappings[mappings[i].read_id].push_back(i);
    }
}

/**
 * Selects the winning node from a window, the winner is decided by the number of matched kmers within the window.
 * A winner is chosen iff the difference between the number of matches for first > second * min_match_spread
 * @param ranking Lookup table of nodes ordered by highest to lowest match count nodes
 * @param min_window_matches Minimum number of matches for winner to be valid
 * @param min_match_spread Minimum spread (match count multiplier) between first and second
 * @return returns the winner's node id, 0 otherwise
 */
int32_t LongReadMapper::getWinner(std::multimap<uint32_t, int32_t> ranking, uint min_window_matches, float min_match_spread) {
    if (ranking.empty()) return 0;
    if (ranking.size() == 1 && ranking.rbegin()->first > min_window_matches) {
        return ranking.rbegin()->second;
    } else if (ranking.size() > 1) {
        auto second=ranking.rbegin();
        ++second;
        if (ranking.rbegin()->first > min_match_spread*second->first ) {
            return ranking.rbegin()->second;
        }
    }
    return 0;
}

uint64_t LongReadMapper::map_reads2(std::string &filename, uint32_t error) {
    FastqReader<FastqRecord> fastqReader({0},filename);
    StrandedMinimiserSketchFactory kf(k, w);
    std::set<MinPosIDX> sketch;
    FastqRecord read;
    while(fastqReader.next_record(read)) {
        std::vector<Match> hits;
        const auto matches = getMatches(read.seq);
        std::copy(matches.begin(),matches.end(), std::ostream_iterator<Match>(sglib::OutputLog(sglib::DEBUG, false)," -> "));
        std::cout << std::endl;

        uint64_t n(matches.size());

        std::vector<int32_t> f(n);
        std::vector<int32_t> t(n);
        std::vector<int32_t> p(n);
        std::vector<int32_t> v(n);
        int64_t st(0);
        uint32_t max_dist_y(5000), max_dist_x(5000);
        uint32_t bw(500);
        uint32_t max_skip(25);
        int64_t min_sc=40;
        int min_cnt=4;
        int32_t i(0), j(0);

        for (i = 0; i < n; ++i){
            uint64_t ri = matches[i].refPos;
            int32_t max_j = -1;
            int32_t qi = matches[i].readPos;
            int32_t n_skip = 0;
            int32_t max_f = this->k;
            int64_t min_d;
            while (st < i && ri - matches[st].refPos > max_dist_x) ++st;
            for (j = i - 1; j >= st; --j) {
                int64_t dr = ri - matches[j].refPos;
                int32_t dq = qi - (int32_t)matches[j].readPos, log_dd;
                int32_t sc;
                int64_t dd;
                if ((dr == 0) || dq <= 0) continue; // don't skip if an anchor is used by multiple segments; see below
                if ((dq > max_dist_y) || dq > max_dist_x) continue;
                dd = dr > dq? dr - dq : dq - dr;
                if (dd > bw) continue;
                if (dr > max_dist_y) continue;
                min_d = dq < dr? dq : dr;
                sc = (int32_t) (min_d > k ? k : dq < dr ? dq : dr);
                log_dd = dd? ilog2_32(dd) : 0;
                sc -= (int)(dd * .01 * k) + (log_dd>>1);
                sc += f[j];
                if (sc > max_f) {
                    max_f = sc, max_j = j;
                    if (n_skip > 0) --n_skip;
                } else if (t[j] == i) {
                    if (++n_skip > max_skip)
                        break;
                }
                if (p[j] >= 0) t[p[j]] = i;
            }
            f[i] = max_f, p[i] = max_j;
            v[i] = max_j >= 0 && v[max_j] > max_f? v[max_j] : max_f; // v[] keeps the peak score up to i; f[] is the score ending at i, not always the peak
        }

        std::fill(t.begin(), t.end(), 0);
        // find the ending positions of chains
        int n_u;
        for (i = 0; i < n; ++i)
            if (p[i] >= 0) t[p[i]] = 1;
        for (i = n_u = 0; i < n; ++i)
            if (t[i] == 0 && v[i] >= min_sc)
                ++n_u;
        if (n_u == 0) {
            continue;
        }

        std::vector<std::pair<int32_t ,int32_t >> u(n_u);
        for (i = n_u = 0; i < matches.size(); ++i) {
            if (t[i] == 0 && v[i] >= min_sc) {
                j = i;
                while (j >= 0 && f[j] < v[j]) j = p[j]; // find the peak that maximizes f[]
                if (j < 0) j = i; // TODO: this should really be assert(j>=0)
                u[n_u] = std::make_pair(f[j],j);
                n_u++;
            }
        }
        std::sort(u.begin(),u.begin()+n_u);
        std::reverse(u.begin(),u.begin()+n_u);

        // backtrack
        std::fill(t.begin(), t.end(), 0);

        int32_t n_v(0), k(0);
        for (i = n_v = k = 0; i < n_u; ++i) { // starting from the highest score
            int32_t n_v0 = n_v, k0 = k;
            j = u[i].second;
            do {
                v[n_v++] = j;
                t[j] = 1;
                j = p[j];
            } while (j >= 0 && t[j] == 0);
            if (j < 0) {
                if (n_v - n_v0 >= min_cnt) {
                    u[k].second = u[i].second;
                    u[k].first = n_v-n_v0;
                    k++;
                }
            } else if (u[i].first - f[j] >= min_sc) {
                if (n_v - n_v0 >= min_cnt) {
                    u[k].second = u[i].second;
                    u[k].first = n_v - n_v0;
                    k++;
                }
            }
            if (k0 == k) n_v = n_v0; // no new chain added, reset
        }

//        *n_u_ = n_u = k, *_u = u; // NB: note that u[] may not be sorted by score here

        // write the result to b[]
        std::vector<Match> b(n_v);
        for (i = 0, k = 0; i < n_u; ++i) {
            int32_t k0 = k, ni = u[i].second;
            for (j = 0; j < ni; ++j)
                b[k] = matches[v[k0 + (ni - j - 1)]], ++k;
        }

        // sort u[] and a[] by a[].x, such that adjacent chains may be joined (required by mm_join_long)
        std::vector<std::pair<uint64_t, uint64_t>> w(n_u);
        for (i = k = 0; i < n_u; ++i) {
            w[i].first = b[k].refPos;
            w[i].second = (uint64_t)k<<32|i;
            k += u[i].second;
        }
        std::sort(w.begin(),w.begin()+n_u);
        std::vector<std::pair<int32_t, int32_t>> u2(n_u * 8);
        std::vector<Match> validMatches;

        for (i = k = 0; i < n_u; ++i) {
            int32_t j = w[i].second, n = u[j].second;
            u2[i] = u[j];
            auto w_i32(w[i].second>>32);
            validMatches.push_back(b[w_i32]);
            k += n;
        }
        std::cout << std::endl;
        std::copy(validMatches.begin(), validMatches.begin()+n_u, std::ostream_iterator<Match> (std::cout, "<==>"));
        std::cout << std::endl;
    }
}

/**
 * Maps a single long read using MinSketches and windowed scoring
 * @param read A single long read in fastq format
 * @return
 * Returns the paths between the nodes in the read, is empty if there are no valid paths across the nodes,
 * if there's no direct link between two nodes, it looks forward until the next node is found. If the next node isn't
 * found after N edges, the path is stored and a new path is started with the next node.
 */
std::vector<LongReadMapping>
LongReadMapper::map_read(FastqRecord read, std::ofstream &matchOutput, std::ofstream &blockOutput) {
    std::vector<LongReadMapping> paths;
    auto matches = getMatchOffsets(read.seq);
    std::sort(matches.begin(),matches.end(), MatchOffset::byReadPos());


    for (const auto &m:matches) {
        matchOutput << "@MATCH " << read.name << ","
                                              << sg.nodes[std::abs(m.dirContig)].sequence.length() << "," << m
                                              << std::endl;
    }

    uint window_size(1000); // 1k window
    uint stepping_size(200);    // 200bp stepping
    auto match_window_multiplier(2.5f);    // MAX ranking node needs to have X times matches in window over second to be winner.
    auto min_window_matches(3u);     // Min times seen a contig in a window to call match
    auto min_windows(5u);    // Min number of windows matching a contig in the same direction
    auto prev = matches.cbegin();
    auto curr = prev;
    std::vector<WindowBlock> blocks;
    for (uint currPos = 0; currPos < read.seq.length(); currPos += stepping_size) {
        std::map<int32_t, uint32_t> nodes;
        if (prev == matches.cend()) break;
        for (curr = prev; curr != matches.cend() and curr->readPos < currPos+window_size; ++curr) {
            nodes[curr->dirContig]++;
        }
        for (; prev != matches.end() and prev->readPos < currPos+stepping_size; ++prev);

        int32_t winner(0);
        if (!nodes.empty()) {
            std::pair<int32_t, uint32_t> max = *nodes.cbegin();
            for (const auto &nc:nodes) {
                if (max.second < nc.second) max = nc;
                sglib::OutputLog(sglib::DEBUG, false) << "@WINDOW," << currPos << "," << currPos + window_size
                                                      << "," << read.name << "," << nc.first << "," << nc.second
                                                      << std::endl;
            }
            std::multimap<uint32_t, int32_t> reverseTest = flip_map(nodes);
            for (std::multimap<uint32_t, int32_t>::const_reverse_iterator it = reverseTest.rbegin();
                 it != reverseTest.rend(); ++it) {
                sglib::OutputLog(sglib::DEBUG, false) << "\t\t@RANK," << currPos << "," << currPos + window_size
                                                      << "," << read.name << ","
                                                      << it->first << "," << it->second << std::endl;
            }

            winner = getWinner(reverseTest, min_window_matches, match_window_multiplier);
        }
        if (!blocks.empty()) blocks.back().count_btw++;
        if (0 == winner) continue;
        if (blocks.empty()) {
            blocks.emplace_back(currPos, currPos+stepping_size, winner, 1);
        }
        else {
            if (blocks.back().contig == winner) {
                blocks.back().end = currPos;
                blocks.back().count++;
            } else {
                blocks.emplace_back(currPos, currPos+stepping_size, winner, 1);
            }
        }
    }

    for (const auto &b: blocks) {
        if (b.count > min_windows) {
            blockOutput << "@BLOCK," << read.name << "," << read.seq.length() << ","
                                    << b.start << "," << b.end << "," << b.contig << "," << b.count
                                    << "," << b.count_btw << std::endl;

            // Find path between blocks, start new SequenceGraphPath if no paths are found.
//            if (!paths.empty()) {
//                LongReadMapping path_between(sg.find_path_between(paths.back().nodes.back(),b.contig, 0, 20));
//                if (!path_between.nodes.empty()) {
//                    for (const auto &node_in_path:path_between.nodes) {
//                        paths.back().append_to_path(node_in_path);
//                    }
//                } else {
//                    paths.push_back(SequenceGraphPath(sg,{b.contig}));
//                }
//            } else {
//                paths.push_back(SequenceGraphPath(sg,{b.contig}));
//            }
        }
    }
    return paths;
}

std::vector<LongReadMapping>
LongReadMapper::map_read(uint32_t readID, std::string sequence, std::ofstream &matchOutput, std::ofstream &blockOutput) {
    std::vector<LongReadMapping> paths;
    auto matches = getMatchOffsets(sequence);
    std::sort(matches.begin(),matches.end(), MatchOffset::byReadPos());


//    for (const auto &m:matches) {
//        matchOutput << "@MATCH " << readID << ","
//                    << sg.nodes[std::abs(m.dirContig)].sequence.length() << "," << m
//                    << std::endl;
//    }

    uint window_size(1000); // 1k window
    uint stepping_size(200);    // 200bp stepping
    auto match_window_multiplier(2.5f);    // MAX ranking node needs to have X times matches in window over second to be winner.
    auto min_window_matches(3u);     // Min times seen a contig in a window to call match
    auto min_windows(5u);    // Min number of windows matching a contig in the same direction
    auto prev = matches.cbegin();
    auto curr = prev;
    std::vector<WindowBlock> blocks;
    for (uint currPos = 0; currPos < sequence.length(); currPos += stepping_size) {
        std::map<int32_t, uint32_t> nodes;
        if (prev == matches.cend()) break;
        for (curr = prev; curr != matches.cend() and curr->readPos < currPos+window_size; ++curr) {
            nodes[curr->dirContig]++;
        }
        for (; prev != matches.end() and prev->readPos < currPos+stepping_size; ++prev);

        int32_t winner(0);
        if (!nodes.empty()) {
            std::pair<int32_t, uint32_t> max = *nodes.cbegin();
            for (const auto &nc:nodes) {
                if (max.second < nc.second) max = nc;
                sglib::OutputLog(sglib::DEBUG, false) << "@WINDOW," << currPos << "," << currPos + window_size
                                                      << "," << readID << "," << nc.first << "," << nc.second
                                                      << std::endl;
            }
            std::multimap<uint32_t, int32_t> reverseTest = flip_map(nodes);
            for (std::multimap<uint32_t, int32_t>::const_reverse_iterator it = reverseTest.rbegin();
                 it != reverseTest.rend(); ++it) {
                sglib::OutputLog(sglib::DEBUG, false) << "\t\t@RANK," << currPos << "," << currPos + window_size
                                                      << "," << readID << ","
                                                      << it->first << "," << it->second << std::endl;
            }

            winner = getWinner(reverseTest, min_window_matches, match_window_multiplier);
        }
        if (!blocks.empty()) blocks.back().count_btw++;
        if (0 == winner) continue;
        if (blocks.empty()) {
            blocks.emplace_back(currPos, currPos+stepping_size, winner, 1);
        }
        else {
            if (blocks.back().contig == winner) {
                blocks.back().end = currPos;
                blocks.back().count++;
            } else {
                blocks.emplace_back(currPos, currPos+stepping_size, winner, 1);
            }
        }
    }

    for (const auto &b: blocks) {
        if (b.count > min_windows) {
            blockOutput << "@BLOCK," << readID << "," << sequence.length() << ","
                        << b.start << "," << b.end << "," << b.contig << "," << b.count
                        << "," << b.count_btw << std::endl;
        }
    }
    return paths;
}

std::vector<LongReadMapper::Match> LongReadMapper::getMatches(std::string &seq) {
    std::vector<Match> matches;
    StrandedMinimiserSketchFactory kf(k, w);
    std::unordered_set<MinPosIDX> sketch;

    kf.getMinSketch(seq, sketch);

    for (const auto &sk:sketch){
        std::unordered_map<uint64_t, std::vector<graphStrandPos>>::const_iterator foundKey(index.find(sk.hash));
        if (foundKey == index.end()) continue;
        for (const auto &match : foundKey->second) {
            matches.emplace_back(match.node * std::signbit(sk.pos)?1:-1, std::abs(sk.pos), match.pos);
        }
    }

    // Sort matches by contig and offset
    std::sort(matches.begin(), matches.end(), typename Match::byContigRef());
    return matches;
}

/**
 * Searches the minSketch database of the reference for each minSketch generated from the query
 * Each match is augmented with the directed node, position on the read and position on the reference
 * @param query The DNA part of the read
 * @return The matching ±node, query position and reference position
 */
std::vector<LongReadMapper::MatchOffset> LongReadMapper::getMatchOffsets(std::string &query) {
    std::vector<MatchOffset> matches;
    uint32_t sketch_not_in_index(0);
    uint32_t sketch_in_index(0);
    StrandedMinimiserSketchFactory kf(k, w);
    std::unordered_set<MinPosIDX> sketch;

    auto read_sketch_num(kf.getMinSketch(query, sketch));

    for (auto sk : sketch) {
        auto foundKey(index.find(sk.hash));
        if (foundKey == index.end()) {
            sketch_not_in_index++;
            continue;
        }
        sketch_in_index++;
        if (foundKey->second.size() < 10) {
            for (auto match : foundKey->second) {
                matches.emplace_back(match.node * (std::signbit(sk.pos) ? -1 : 1), std::abs(sk.pos), match.pos);
            }
        }
    }

    std::sort(matches.begin(), matches.end(), typename MatchOffset::byOffset());
    sglib::OutputLog(sglib::DEBUG, false) << "Matched " << sketch_in_index
                                          << " / " << read_sketch_num << " read sketches" << std::endl;
    return matches;
}
