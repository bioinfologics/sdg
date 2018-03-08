//
// Created by Luis Yanes (EI) on 15/09/2017.
//

#ifndef SG_CONTIGBLOCK_H_H
#define SG_CONTIGBLOCK_H_H

#include <fstream>
#include <iostream>
#include <unordered_set>
#include <set>
#include <limits>
#include <iomanip>
#include <cmath>
#include <numeric>
#include <unordered_map>
#include <algorithm>
#include <iterator>
#include <map>
#include <sglib/logger/OutputLog.h>
#include <sglib/readers/Common.h>
#include "KMerIDXFactory.h"

struct FilterSetParams {
    FilterSetParams(std::string output_prefix, uint8_t k, std::unordered_set<KmerIDX> &uniq_kmers,
                    uint32_t min_kmers_to_call_match = 1,
                    uint32_t min_seen_contig_to_write_output = 0) : k(k), uniq_kmers(uniq_kmers),
                                                                    min_kmers_to_call_match(min_kmers_to_call_match),
                                                                    output_prefix(output_prefix),
                                                                    min_seen_contig_to_write_output(
                                                                            min_seen_contig_to_write_output) {}

    uint8_t k;
    const std::unordered_set<KmerIDX> &uniq_kmers;
    const uint32_t min_kmers_to_call_match;
    const uint32_t min_seen_contig_to_write_output;
    const std::string output_prefix;
    const std::unordered_map<int32_t, std::pair<uint64_t, uint32_t>> unitigs;
};

template <typename FileRecord>
class ContigBlockFactory : protected KMerFactory {
    template<typename T, typename T2=T>
    struct accumulator {
        T2 sum; // we could plug in a more accurate type for the sum
        T S;
        T M;
        size_t N;

        // default constructor initializes all values
        accumulator() : sum(0), S(0), M(0), N(0) {}

        // add another number
        T2 operator()(const T &x) {
            ++N;
            sum += x;
            T Mprev = M;
            M += (x - Mprev) / N;
            S += (x - Mprev) * (x - M);
            return sum;
        }

        T mean() const { return sum / N; }

        T variance() const { return S / (N - 1); }

        // operator<< to print the statistics to screen:
        // denoted friend just to be able to write this inside
        // the class definition and thus not to need to write
        // the template specification of accumulator...
        friend std::ostream &operator<<(std::ostream &out, const accumulator &a) {
            if (a.N > 0)
                out << "N\t\t\t= " << a.N << std::endl << "sum\t\t\t= " << a.sum << std::endl << "mean\t\t= "
                    << std::fixed << std::setprecision(2) << a.mean() << std::endl;
            if (a.N > 1)
                out << "sd\t\t\t= " << std::fixed << std::setprecision(2) << std::sqrt(a.variance()) << std::endl;
            else
                out << "sd\t\t\t= " << std::fixed << std::setprecision(2) << 0 << std::endl;
            return out;
        }

    };

public:
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

        friend std::ostream &operator<<(std::ostream &os, const Block &block) {
            os << block.contigID << " ref("
               << block.rstart << ";" << block.rend << ") query("
               << block.qstart << ";" << block.qend << "):"
               << block.count;
            return os;
        }

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

    };
    struct Match {
    public:
        Match() : dirContig(0), readPos(0), refPos(0) {};

        Match(int32_t contig, uint32_t readPos, uint32_t refPos) :
                dirContig(contig), readPos(readPos), refPos(refPos) {}

        int32_t dirContig;  // Sign indicates direction
        uint32_t readPos;   // Pos on the read
        uint32_t refPos;    // Ref pos

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

    const bool isValid(const Block &b) const {

        //  Minimo numero de kmers en el bloque
        if (b.count < min_kmers_to_call_match) {
            return false;
        }
        return true;
    }

    explicit ContigBlockFactory(FilterSetParams params) :
            KMerFactory(params.k),
            bases(0),
            kmers(params.uniq_kmers),
            min_kmers_to_call_match(params.min_kmers_to_call_match),
            min_seen_contig_to_write_output(params.min_seen_contig_to_write_output), offset_limit(1000),
            ctg_read(params.output_prefix+"ctg_read.txt"),
            read_contig_link(params.output_prefix + "read_contig_link.txt"),
            read_length(params.output_prefix + "read_length.log"),
            read_link_stats(params.output_prefix + "read_link.log"),
            read_block_stats(params.output_prefix + "read_block.log"),
            read_coverage_stats(params.output_prefix + "read_coverage.log")
    {
        read_block_stats << "blockCtgDif, blockOffsetDif, Validblock, Tblock" << std::endl;
        read_coverage_stats << "Mappedkmers, validBlockKmers, Tkmers" << std::endl;
        read_length << "Name, Length" << std::endl;
        read_link_stats << "Number of links" << std::endl;
        read_contig_link << "Links" << std::endl;
    };
    ~ContigBlockFactory() = default;

    void setFileRecord(FileRecord& rec){
        currentFileRecord = rec;
        fkmer=0;
        rkmer=0;
        last_unknown=0;
    }


    std::vector<Block> getBlocks(const std::vector<Match> &matches) const {
        std::vector<Block> blocks;
        uint32_t blkDif(0), offDif(0);
        if (!matches.empty()) {
            auto prev = matches.cbegin();
            Block blk(0, prev->refPos, 0, prev->readPos, 0, 1, prev->dirContig);
            auto curr = prev;
            ++curr;
            for (; curr != matches.cend(); ++curr) {
                int32_t currOffset(curr->refPos-curr->readPos);
                int32_t prevOffset(prev->refPos - prev->readPos);
                if (blk.contigID == curr->dirContig
                    and std::abs(currOffset - prevOffset) < offset_limit ) {
                    blk.count++;
                    blk.setREnd(curr->refPos), blk.setQEnd(curr->readPos);
                } else {
                    if (blk.contigID != curr->dirContig) blkDif++;
                    if (std::abs((int32_t)(curr->refPos - prev->refPos)) >= offset_limit) offDif++;
                    if (blk.count>min_kmers_to_call_match)
                        blocks.push_back(blk);
                    blk = Block(curr - matches.begin(), curr->refPos, 0, curr->readPos, 0, 1, curr->dirContig);
                }
                prev = curr;
            }
            if(blk.count > min_kmers_to_call_match)
                blocks.push_back(blk);
        }
        read_block_stats << blkDif << "," << offDif; // First 2 values, missing validBlocks
        return blocks;
    }

    std::vector<Match>  getMatches() {
        std::vector<Match> matches;
        int64_t offset(0);
        uint64_t p(0);
        while (p < currentFileRecord.seq.size()) {
            //fkmer: grows from the right (LSB)
            //rkmer: grows from the left (MSB)
            fillKBuf(currentFileRecord.seq[p], fkmer, rkmer, last_unknown);
            p++;
            if (last_unknown >= K) {
                // Generate all the kmer-contig matches with Direction, Contig, Offset and readPosition
                std::pair<KmerIDX,int> read_kmer (
                        (fkmer < rkmer) ? std::make_pair(KmerIDX(fkmer), 1) : std::make_pair(KmerIDX(rkmer), -1));
                auto kmer = kmers.find(read_kmer.first);
                if (kmer != kmers.end()) { // If the kmer from the read is in the unitigs
                    if (kmer->contigID * read_kmer.second > 0) {
                        offset = kmer->pos - p;
                    } else {
                        offset = kmer->pos + p;
                    }
                    matches.emplace_back(kmer->contigID * read_kmer.second, p, kmer->pos);
                }
                ctg_read << ((kmer != kmers.end()) ? abs(kmer->contigID) : 0) << " ";
            }
        }
        // Sort matches by contig and offset
        std::sort(matches.begin(), matches.end(), typename Match::byContigRef());
        return matches;
    }

    const std::vector<std::vector<std::pair<int32_t, unsigned int>>> getTop(unsigned int n, std::vector<Match> &matches, unsigned int read_length, unsigned int window) const {
        std::vector<std::vector<std::pair<int32_t, unsigned int>>> result(1+ read_length/window);
        if (matches.empty()) return result;
        auto pos (matches[0].readPos);

        typename std::vector<Match>::const_iterator m(matches.cbegin());
        std::vector<std::map<int32_t, unsigned int>> contigMatch_section(1 + read_length/window);
        for (; m != matches.cend(); ++m){
            contigMatch_section[m->readPos/window][m->dirContig]++;
        }

        for (auto s = 0; s!=read_length/window; s++) {
            std::vector<std::pair<int32_t , unsigned int>> top(n);
            std::partial_sort_copy(contigMatch_section[s].begin(),
                                   contigMatch_section[s].end(),
                                   top.begin(),
                                   top.end(),
                                   [](std::pair<const int32_t, unsigned int> const& l,
                                      std::pair<const int32_t, unsigned int> const& r)
                                   {
                                       return l.second > r.second;
                                   });
            std::cout << "Section " << s << "\t";
            for (const auto &t:top) {
                std::cout << t.first << ':' << t.second << " | ";
            }
            std::cout << std::endl;
            std::copy(top.begin(),top.end(), std::back_inserter(result[s]));
        }
        return result;
    }


    const bool next_element(std::vector<Block> &blocks) {
        blocks.clear();
        std::cout << currentFileRecord.name << "\t" << currentFileRecord.seq.size() << std::endl;
        bases += currentFileRecord.seq.size();
        read_length << currentFileRecord.name << "," << currentFileRecord.seq.size() << std::endl;
        ctg_read << currentFileRecord.name << "(" << currentFileRecord.seq.size() << ")\n";
        std::vector<Match> matches = getMatches();  // Return all matches sorted by position in the read
        ctg_read << std::endl;
        std::copy(matches.begin(), matches.end(),
                  std::ostream_iterator<Match>(sglib::OutputLog(sglib::DEBUG, false), " -> "));
        std::cout << std::endl;

        // Generate blocks which share the same ctgOffset
        blocks = getBlocks(matches);
        std::cout << std::endl;
        std::sort(blocks.begin(),blocks.end(), typename Block::byRefPos());
//         Sort blocks by ReadPos and keep only the valid ones
        std::copy(blocks.cbegin(), blocks.cend(), std::ostream_iterator<Block>(std::cout, ", "));
        std::cout << std::endl;

        // print the rest of the parameters for the blocks on this read
        read_block_stats << "," << blocks.size() << "," << blocks.size() << std::endl;
        std::cout << std::endl;
        read_link_stats << blocks.size() << std::endl;
        return false;
    }


    const bool next_elemen2t(std::vector<Block> &blocks) {
        blocks.clear();
        std::cout << currentFileRecord.name << "\t" << currentFileRecord.seq.size() << std::endl;
        bases+=currentFileRecord.seq.size();
        read_length << currentFileRecord.name << "," << currentFileRecord.seq.size() << std::endl;
        ctg_read << currentFileRecord.name << "(" << currentFileRecord.seq.size() << ")\n";
        std::vector<Match> matches = getMatches();  // Return all matches sorted by position in the read
        ctg_read<<std::endl;
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
            int32_t max_f = this->K;
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
                sc = (int32_t) (min_d > K ? K : dq < dr ? dq : dr);
                log_dd = dd? ilog2_32(dd) : 0;
                sc -= (int)(dd * .01 * K) + (log_dd>>1);
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
            return false;
        }

        std::vector<uint64_t > u(static_cast<unsigned long>(n_u * 8));
        for (i = n_u = 0; i < matches.size(); ++i) {
            if (t[i] == 0 && v[i] >= min_sc) {
                j = i;
                while (j >= 0 && f[j] < v[j]) j = p[j]; // find the peak that maximizes f[]
                if (j < 0) j = i; // TODO: this should really be assert(j>=0)
                u[n_u++] = (uint64_t)f[j] << 32 | j;
            }
        }
        std::sort(u.begin(),u.begin()+n_u);
        for (i = 0; i < n_u>>1; ++i) { // reverse, s.t. the highest scoring chain is the first
            uint64_t t = u[i];
            u[i] = u[n_u - i - 1], u[n_u - i - 1] = t;
        }

        // backtrack
        std::fill(t.begin(), t.end(), 0);

        int32_t n_v(0), k(0);
        for (i = n_v = k = 0; i < n_u; ++i) { // starting from the highest score
            int32_t n_v0 = n_v, k0 = k;
            j = (int32_t)u[i];
            do {
                v[n_v++] = j;
                t[j] = 1;
                j = p[j];
            } while (j >= 0 && t[j] == 0);
            if (j < 0) {
                if (n_v - n_v0 >= min_cnt) u[k++] = u[i]>>32<<32 | (n_v - n_v0);
            } else if ((int32_t)(u[i]>>32) - f[j] >= min_sc) {
                if (n_v - n_v0 >= min_cnt) u[k++] = ((u[i]>>32) - f[j]) << 32 | (n_v - n_v0);
            }
            if (k0 == k) n_v = n_v0; // no new chain added, reset
        }

//        *n_u_ = n_u = k, *_u = u; // NB: note that u[] may not be sorted by score here

        // write the result to b[]
        std::vector<Match> b(n_v);
        for (i = 0, k = 0; i < n_u; ++i) {
            int32_t k0 = k, ni = (int32_t)u[i];
            for (j = 0; j < ni; ++j)
                b[k] = matches[v[k0 + (ni - j - 1)]], ++k;
        }

        // sort u[] and a[] by a[].x, such that adjacent chains may be joined (required by mm_join_long)
        std::vector<std::pair<uint64_t, uint64_t>> w(n_u);
        for (i = k = 0; i < n_u; ++i) {
            w[i].first = b[k].refPos, w[i].second = (uint64_t)k<<32|i;
            k += (int32_t)u[i];
        }
        std::sort(w.begin(),w.begin()+n_u);
        std::vector<uint64_t > u2(n_u * 8);
        std::vector<Match> validMatches;
        for (i = k = 0; i < n_u; ++i) {
            int32_t j = (int32_t)w[i].second, n = (int32_t)u[j];
            u2[i] = u[j];
            validMatches.push_back(b[w[i].second>>32]);
            k += n;
        }

        // Get contigs by "appearances" and "kmer matches" within 1kbp window for every 1kbp window of the read

        auto tmpTop = getTop(5, matches, currentFileRecord.seq.length(), 1000);

//        std::sort(matches.begin(),matches.rend(), typename Match::byReadPos());
        // Generate blocks which share the same ctgOffset
        auto tmpBlocks = getBlocks(matches);
        std::copy(tmpBlocks.cbegin(), tmpBlocks.cend(), std::ostream_iterator<Block>(std::cout, ", "));
        std::cout << std::endl;

//         Sort blocks by ReadPos and keep only the valid ones
        std::copy_if(tmpBlocks.begin(), tmpBlocks.end(), std::back_inserter(blocks),
                     [&](const Block &b) { return isValid(b); });
        std::copy(blocks.cbegin(), blocks.cend(), std::ostream_iterator<Block>(std::cout, ", "));
        std::cout << std::endl;


        // print the rest of the parameters for the blocks on this read
        read_block_stats << "," << blocks.size() << "," << blocks.size() << std::endl;
        std::cout << std::endl;
        read_link_stats << blocks.size() << std::endl;
        return false;
    }

private:

    uint64_t bases;

    uint32_t currentContigID;
    FileRecord currentFileRecord;

    std::ofstream read_contig_link; /// INFO - Links seen so far

    const uint32_t min_kmers_to_call_match;
    const uint32_t min_seen_contig_to_write_output;
    const uint32_t offset_limit;


    const std::unordered_set<KmerIDX>& kmers;

    char b2f[255];
    char b2r[255];

    // T = Total
    std::ofstream ctg_read;                     /// DEBUG: output only for debug mode
    mutable std::ofstream read_length;          /// INFO - Read, Length
    mutable std::ofstream read_coverage_stats;  /// INFO - Mapped kmers, ValidBlockKmers, TKmer
    mutable std::ofstream read_block_stats;     /// INFO - blockCtgDif, blockOffsetDif, Validblock, Tblock
    mutable std::ofstream read_link_stats;      /// INFO - # Links generated
};

#endif //SG_CONTIGBLOCK_H_H
