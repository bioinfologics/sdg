//
// Created by Luis Yanes (EI) on 12/02/2018.
//

#ifndef BSG_LONGREADMAPPER_H
#define BSG_LONGREADMAPPER_H


#include <iostream>
#include <sglib/factories/KMerIDXFactory.h>
#include <sglib/readers/FileReader.h>
#include <sglib/readers/SequenceGraphReader.h>
#include <sglib/SMR.h>
#include <sglib/PairedReadMapper.h>
#include <sglib/factories/StrandedMinSketchFactory.h>
#include <sglib/indexers/minSketchIndex.hpp>

template <typename A, typename B>
std::multimap<B, A> flip_map(std::map<A,B> & src) {

    std::multimap<B,A> dst;

    for(typename std::map<A, B>::const_iterator it = src.begin(); it != src.end(); ++it)
        dst.insert(std::pair<B, A>(it -> second, it -> first));

    return dst;
}

class LongReadMapper {

    SequenceGraph & sg;
    minSketchIndex index;
    uint8_t k=15;
    uint8_t w=5;
    std::unordered_map<uint64_t, std::vector<graphStrandPos>> kmer_to_graphposition;

    struct MatchOffset{
        MatchOffset(){};
        MatchOffset(int32_t dirContig, uint32_t readPos, uint32_t offset) :
                dirContig(dirContig), readPos(readPos), offset(offset) {};

        friend std::ostream &operator<<(std::ostream &os, const MatchOffset &match) {
            os << "" << match.dirContig << ","
               << match.readPos << ","
               << match.offset
               << "";
            return os;
        }
        int32_t dirContig;  // Sign indicates direction
        int32_t readPos;   // Pos on the read
        int32_t offset;    // Ref offset (read_start vs ref_start)

        friend class byOffset;
        struct byOffset {
            bool operator()(const MatchOffset &a, const MatchOffset &b) {
                return std::tie(a.offset) < std::tie(b.offset);
            }
        };

        friend class byContigRead;
        struct byContigRead {
            bool operator()(const MatchOffset &a, const MatchOffset &b) {
                return std::tie(a.dirContig, a.readPos) < std::tie(b.dirContig, b.readPos);
            }
        };

        friend class byContigOffsetPos;
        struct byContigOffsetPos {
            bool operator()(const MatchOffset &a, const MatchOffset &b) {
                return std::tie(a.dirContig, a.offset, a.readPos) < std::tie(b.dirContig, b.offset, b.readPos);
            }
        };

        friend class byReadPos;
        struct byReadPos {
            bool operator()(const MatchOffset &a, const MatchOffset &b) {
                return std::tie(a.readPos) < std::tie(b.readPos);
            }
        };
    };
    struct Match {
    public:
        Match() : dirContig(0), readPos(0), refPos(0) {};

        Match(int32_t contig, uint32_t readPos, uint32_t refPos) :
                dirContig(contig), readPos(readPos), refPos(refPos) {}

        int32_t dirContig;  // Sign indicates direction
        int32_t readPos;   // Pos on the read
        int32_t refPos;    // Ref pos

        friend std::ostream &operator<<(std::ostream &os, const Match &match) {
            os << "(" << match.dirContig << "," << match.refPos
               //               << ","
               //               << match.readPos
               << ")";
            return os;
        }


        friend class byReadPos;
        struct byReadPos {
            bool operator()(const Match &a, const Match &b) {
                return std::tie(a.readPos) < std::tie(b.readPos);
            }
        };

        friend class byContigRef;
        struct byContigRef {
            bool operator()(const Match &a, const Match &b) {
                return std::tie(a.dirContig,a.refPos) < std::tie(b.dirContig,b.refPos);
            }
        };
    };
    struct Block {
        Block() : rstart(0), rend(0), contigID(0), count(0) {};

        Block(size_t readID, uint32_t start, uint32_t end, uint32_t qstart, uint32_t qend, uint16_t count,
              int32_t contig) : readID(readID),
                                rstart(start), rend(end),
                                qstart(qstart), qend(qend),
                                contigID(contig),
                                count(count) {}

        void setREnd(uint32_t val) {
            rend = val;
        }

        void setQEnd(uint32_t val) {
            qend = val;
        }

        size_t readID = 0;
        uint32_t rstart;     // Block match ref_start
        uint32_t rend;       // Block match ref_end
        uint32_t qstart;    // Query start
        uint32_t qend;      // Query end
        int32_t contigID;   // Sign indicates direction
        uint16_t count;     // Kmer matches in the block


        friend class byCount;
        struct byCount {
            bool operator()(const Block &a, const Block &b) {
                return std::tie(a.count,a.rstart) > std::tie(b.count,b.rstart);
            }
        };

        friend class byRefPos;
        struct byRefPos {
            bool operator()(const Block &a, const Block &b) {
                return std::tie(a.rstart) < std::tie(b.rstart);
            }
        };

        friend class byQuery;
        struct byQuery {
            bool operator()(const Block &a, const Block &o) {
                return std::tie(a.contigID, a.rstart ) < std::tie(o.contigID, o.rstart);
            }
        };

        friend class byContig;
        struct byContig {
            bool operator()(const Block &a, const Block &o) {
                return std::tie(a.contigID) < std::tie(o.contigID);
            }
        };

        bool operator==(const Block &o) const {return std::tie(readID, contigID) == std::tie(o.readID,o.contigID);}
        bool operator<(const Block &o) const {return std::tie(readID, contigID,rstart,rend) < std::tie(o.readID,o.contigID,rstart,rend);}
        bool operator>(const Block &o) const {return std::tie(readID, contigID,rstart,rend) < std::tie(o.readID,o.contigID,rstart,rend);}

        void merge(const Block &o) {
            if (rstart > o.rstart) {rstart = o.rstart;}
            if (rend < o.rend){rend = o.rend;}
            count++;
        }

        friend std::istream& operator>>(std::istream& is, const Block& block) {
            is.read((char*)&block, sizeof(block));
            return is;
        }

        friend std::ostream& operator<<(std::ostream& os, const Block& block) {
            os << block.qstart << "," << block.qend << "," << (std::signbit(block.contigID)?'-':'+' ) << "," << std::abs(block.contigID) <<"," << block.rstart << "," << block.rend;
            return os;
        }
    };
    struct WindowBlock {
        WindowBlock(uint32_t s, uint32_t e, int32_t contig, uint32_t count) : start(s), end(e), contig(contig), count(count) {}
        uint32_t start = 0;
        uint32_t end = 0;
        int32_t contig = 0;
        uint32_t count = 1;
        uint32_t count_btw = 1;
    };

    std::vector<Match> getMatches(std::string &seq) {
        std::vector<Match> matches;
        int64_t offset(0);
        uint64_t p(0);
        StrandedMinimiserSketchFactory kf(k, w);
        std::set<MinPosIDX> sketch;

        kf.getMinSketch(seq, sketch);

        for (const auto &sk:sketch){
            std::unordered_map<uint64_t, std::vector<graphStrandPos>>::const_iterator foundKey(kmer_to_graphposition.find(sk.hash));
            if (foundKey == kmer_to_graphposition.end()) continue;
            for (const auto &match : foundKey->second) {
                matches.emplace_back(match.node * std::signbit(sk.pos)?1:-1, std::abs(sk.pos), std::abs(match.pos));
            }
        }

        // Sort matches by contig and offset
        std::sort(matches.begin(), matches.end(), typename Match::byContigRef());
        return matches;
    }

    std::vector<MatchOffset> getMatchOffsets(std::string &seq) {
        std::vector<MatchOffset> matches;
        uint32_t sketch_not_in_index(0);
        uint32_t sketch_in_index(0);
        int64_t offset(0);
        uint64_t p(0);
        StrandedMinimiserSketchFactory kf(k, w);
        std::set<MinPosIDX> sketch;

        auto read_sketch_num(kf.getMinSketch(seq, sketch));

        for (auto sk : sketch) {
            auto foundKey(index.find(sk.hash));
            if (foundKey == index.end()) {
                sketch_not_in_index++;
                continue;
            }
            sketch_in_index++;
            for (auto match : foundKey->second) {
                matches.emplace_back(match.node*(std::signbit(sk.pos)?-1:1), std::abs(sk.pos), std::abs(match.pos));
            }
        }

        std::sort(matches.begin(), matches.end(), typename MatchOffset::byOffset());
        sglib::OutputLog(sglib::DEBUG, false) << "Matched " << sketch_in_index
                                             << " / " << read_sketch_num << " read sketches" << std::endl;
        return matches;
    }

public:
    LongReadMapper(uint8_t k, uint8_t w, SequenceGraph &sg) : sg(sg), k(k), w(w), index(sg, k, w) {}

    void map_read(FastqRecord read) {
        auto matches = getMatchOffsets(read.seq);
        std::sort(matches.begin(),matches.end(), MatchOffset::byReadPos());

        for (const auto &m:matches) {
            sglib::OutputLog(sglib::DEBUG, false) << "@MATCH " << read.name << ","
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
                sglib::OutputLog(false) << "@BLOCK," << read.name << "," << read.seq.length() << ","
                                                     << b.start << "," << b.end << "," << b.contig << "," << b.count
                                                     << "," << b.count_btw << std::endl;
            }
        }
    }

    uint64_t map_reads2(std::string &filename, uint32_t error) {
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

    int32_t getWinner(std::multimap<uint32_t , int32_t> ranking, uint min_window_matches, float min_match_spread) {
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

    uint64_t map_reads(std::string &filename, uint32_t error) {

        FastqReader<FastqRecord> fastqReader({0}, filename);
        StrandedMinimiserSketchFactory kf(k, w);
        std::set<MinPosIDX> sketch;
        FastqRecord read;

        std::cout << "@MATCH,READ,NODE_LENGTH,NODE,READ_POS,OFFSET" << std::endl;
        std::cout << "@BLOCK,READ,READ_LENGTH,START,END,NODE,COUNT,COUNT_BTW" << std::endl;
        while (fastqReader.next_record(read)) {
            map_read(read);
        }
        return 0;
    }
};


#endif //BSG_LONGREADMAPPER_H
