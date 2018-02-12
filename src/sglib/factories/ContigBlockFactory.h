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

struct Block {
    Block() : start(0), end(0), contigID(0), offset(0), count(0) {};

    Block(size_t readID, uint32_t start, int32_t contig, int64_t offset, uint16_t count, uint32_t end) : readID(readID), start(start), end(end),
                                                                                          contigID(contig),
                                                                                          offset(offset),
                                                                                          count(count) {
        if (contigID < 0) {
            this->start = end;
            this->end = start;
        }
    }

    void setEnd(uint32_t val) {
        if (contigID < 0)
            start = val;
        else
            end = val;
    }

    size_t readID;
    uint32_t start;     // Block match start
    uint32_t end;       // Block match end
    int32_t contigID;   // Sign indicates direction
    int64_t offset;     // Pos_contig - Pos_read
    uint16_t count;     // Kmer matches in the block

    friend std::ostream &operator<<(std::ostream &os, const Block &block) {
        os << block.contigID << " (" << block.start << ";" << block.end << ";" << int(block.end - block.start)
           << ";" << block.count << ")";
        return os;
    }

    friend class byReadPos;
    friend class byCount;
    struct byCount {
        bool operator()(const Block &a, const Block &b) {
            return a.count < b.count;
        }
    };
    struct byReadPos {
        bool operator()(const Block &a, const Block &o) {
            return std::tie(a.start, a.contigID) < std::tie(o.start, o.contigID);
        }
    };

    bool operator==(const Block &o) const {return std::tie(readID, contigID) == std::tie(o.readID,o.contigID);}
    bool operator<(const Block &o) const {return std::tie(readID, contigID,start,end) < std::tie(o.readID,o.contigID,start,end);}
    bool operator>(const Block &o) const {return std::tie(readID, contigID,start,end) < std::tie(o.readID,o.contigID,start,end);}

    void merge(const Block &o) {
        if (start > o.start) {start = o.start;}
        if (end < o.end){end = o.end;}
        if (std::abs(offset) > std::abs(o.offset)) offset = o.offset;
        count++;
    }

    friend std::istream& operator>>(std::istream& is, const Block& block) {
        is.read((char*)&block, sizeof(block));
        return is;
    }

};

struct Match {
public:
    Match() : dirContig(0), offset(0), readPos(0) {};

    Match(int32_t contig, int64_t offset, uint32_t readPos) : dirContig(contig), offset(offset), readPos(readPos) {}

    int32_t dirContig;  // Sign indicates direction
    int64_t offset;     // Pos_contig - Pos_read
    uint32_t readPos;

    friend std::ostream &operator<<(std::ostream &os, const Match &match) {
        os << "(" << match.dirContig << "," << match.offset
           //               << ","
           //               << match.readPos
           << ")";
        return os;
    }

    friend class byCtgOffset;

    struct byCtgOffset {
        bool operator()(const Match &a, const Match &b) {
            return std::tie(a.dirContig, a.offset) < std::tie(b.dirContig, b.offset);
        }
    };
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
            Block blk(0, prev->readPos, prev->dirContig, prev->offset, 1, 0);
            auto curr = prev;
            ++curr;
            for (; curr != matches.cend(); ++curr) {
                if (blk.contigID == curr->dirContig and abs(curr->offset - prev->offset) < offset_limit) {
                    blk.count++;
                    blk.setEnd(curr->readPos);
                } else {
                    if (blk.contigID != curr->dirContig) blkDif++;
                    if (abs(curr->offset - prev->offset) >= offset_limit) offDif++;
                    blocks.push_back(blk);
                    blk = Block(matches.begin()-curr, curr->readPos, curr->dirContig, curr->offset, 1, 0);
                }
                prev = curr;
            }
            blocks.push_back(blk);
        }
        std::sort(blocks.begin(), blocks.end(), typename Block::byReadPos());
//        read_block_stats << blkDif << "," << offDif; // First 2 values, missing validBlocks
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
                        offset = kmer->pos - currentFileRecord.seq.size() + p;
                    }
                    matches.emplace_back(kmer->contigID * read_kmer.second, offset, p);
                }
//                ctg_read << ((kmer != kmers.end()) ? abs(kmer->contigID) : 0) << " ";
            }
        }
        // Sort matches by contig and offset
        std::sort(matches.begin(), matches.end(), typename Match::byCtgOffset());
        return matches;
    }

    const bool next_element(std::vector<Block> &blocks) {
        blocks.clear();
//        std::cout << currentFileRecord.name << "\t" << currentFileRecord.seq.size() << " ";
        bases+=currentFileRecord.seq.size();
//        read_length << currentFileRecord.name << "," << currentFileRecord.seq.size() << std::endl;
//        ctg_read << ">" << currentFileRecord.name << "(" << currentFileRecord.seq.size() << ")\n";
//        ctg_read << std::endl;

        std::vector<Match> matches = getMatches();
//        std::copy(matches.begin(),matches.end(),std::ostream_iterator<Match>(std::cout, ";"));std::cout << std::endl;

        // Generate blocks which share the same ctgOffset
        blocks = getBlocks(matches);

        // Sort blocks by ReadPos and keep only the valid ones
//        std::copy_if(blocks.begin(), blocks.end(), std::back_inserter(sortedMatchbyCtgOffset),
//                     [&](const Match &b) { return isValid(b); });
//        std::copy(validBlocks.cbegin(), validBlocks.cend(), std::ostream_iterator<Block>(std::cout, ";"));

        // print the rest of the parameters for the blocks on this read
//        read_block_stats << "," << validBlocks.size() << "," << blocks.size() << std::endl;
//        std::cout << std::endl;
//        read_link_stats << blocks.size() << std::endl;
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
