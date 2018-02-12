//
// Created by Luis Yanes (EI) on 15/09/2017.
//

#ifndef SEQSORTER_CONTIGLINK_H_H
#define SEQSORTER_CONTIGLINK_H_H

#include <fstream>
#include <iostream>
#include <unordered_set>
#include <set>
#include <limits>
#include <iomanip>
#include <cmath>
#include <numeric>
#include "KMerIDXFactory.h"

struct FilterSetParams {
    FilterSetParams(std::string output_prefix, uint8_t k, std::vector<KmerIDX> &uniq_kmers,
                    std::unordered_map<int32_t, std::pair<uint64_t, uint32_t>> &sequences, uint32_t min_kmers_to_call_match = 1,
                    uint32_t min_seen_contig_to_write_output = 0) : k(k), uniq_kmers(uniq_kmers),
                                                                    min_kmers_to_call_match(min_kmers_to_call_match),
                                                                    unitigs(sequences), output_prefix(output_prefix),
                                                                    min_seen_contig_to_write_output(
                                                                            min_seen_contig_to_write_output) {}

    uint8_t k;
    const std::vector<KmerIDX> uniq_kmers;
    const uint32_t min_kmers_to_call_match;
    const uint32_t min_seen_contig_to_write_output;
    const std::string output_prefix;
    const std::unordered_map<int32_t, std::pair<uint64_t, uint32_t>> unitigs;
};

class ContigLink {
public:
    ContigLink(): contig1(std::numeric_limits<int32_t>::max()), contig2(std::numeric_limits<int32_t>::max()), count(std::numeric_limits<uint16_t>::max()) {}

    ContigLink(int32_t c1, int32_t c2, uint16_t _count = 1) : contig1(c1), contig2(c2), count(_count) {}
    int32_t contig1;
    int32_t contig2;
    uint16_t count;

    const bool operator<(const ContigLink& other) const {
        if (contig1<other.contig1) return true;
        if (contig1>other.contig1) return false;
        return contig2<other.contig2;
    }

    const bool operator>(const ContigLink& other) const {
        return ( (contig1 > other.contig1) and (contig2 > other.contig2) );
    }

    const bool operator==(const ContigLink &other) const {
        return ( (contig1 == other.contig1) and (contig2 == other.contig2) );
    }
    void merge(const ContigLink &other) {
        count += other.count;
    }

    ContigLink max() { return ContigLink(); }
    friend std::ostream& operator<<(std::ostream& os, const ContigLink& link){
        os << std::setw(5) << link.contig1 << "<>" << std::setw(5) << link.contig2 << ":" << (int) link.count;
        return os;
    }

    struct byCount {
        bool operator()(const ContigLink &a, const ContigLink &b) {
            return (a.count < b.count);
        }
    };

    friend std::istream& operator>>(std::istream& is, const ContigLink& link) {
        is.read((char*)&link, sizeof(link));
        return is;
    }

};

template <typename FileRecord>
class ContigLinkFactory : protected KMerFactory {
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
    struct Block {
        Block() : start(0), end(0), contigID(0), offset(0), count(0) {};

        Block(uint32_t start, int32_t contig, int64_t offset, uint16_t count, uint32_t end) : start(start), end(end),
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

        struct byReadPos {
            bool operator()(const Block &a, const Block &o) {
                return std::tie(a.start, a.contigID) < std::tie(o.start, o.contigID);
            }
        };
    };

    const bool isValid(const Block &b) const {
        const auto seq_contig = unitigs.find(std::abs(b.contigID));
        if (seq_contig == unitigs.end()) return false;

        //  Minimo numero de kmers en el bloque
        if (b.count < min_kmers_to_call_match) {
            return false;
        }
        return true;
    }

    explicit ContigLinkFactory(FilterSetParams params) :
            KMerFactory(params.k),
            bases(0),
            kmers(params.uniq_kmers.cbegin(), params.uniq_kmers.cend()),
            unitigs(params.unitigs),
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
    ~ContigLinkFactory() = default;

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
            Block blk(prev->readPos, prev->dirContig, prev->offset, 1, 0);
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
                    blk = Block(curr->readPos, curr->dirContig, curr->offset, 1, 0);
                }
                prev = curr;
            }
            blocks.push_back(blk);
        }
        std::sort(blocks.begin(), blocks.end(), typename Block::byReadPos());
        read_block_stats << blkDif << "," << offDif; // First 2 values, missing validBlocks
        return blocks;
    }

    std::vector<Match> getMatches() {
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

    const bool next_element(std::vector<ContigLink> &links) {
        links.clear();
        std::cout << currentFileRecord.name << "\t" << currentFileRecord.seq.size() << " ";
        bases+=currentFileRecord.seq.size();
        read_length << currentFileRecord.name << "," << currentFileRecord.seq.size() << std::endl;
//        ctg_read << ">" << currentFileRecord.name << "(" << currentFileRecord.seq.size() << ")\n";
//        ctg_read << std::endl;

        std::vector<Match> matches = getMatches();
//        std::copy(matches.begin(),matches.end(),std::ostream_iterator<Match>(std::cout, ";"));std::cout << std::endl;

        // Generate blocks which share the same ctgOffset
        std::vector<Block> blocks = getBlocks(matches);

        // Sort blocks by ReadPos and keep only the valid ones
        std::vector<Block> validBlocks;

        std::copy_if(blocks.begin(), blocks.end(), std::back_inserter(validBlocks),
                     [&](const Block &b) { return isValid(b); });
        std::copy(validBlocks.cbegin(), validBlocks.cend(), std::ostream_iterator<Block>(std::cout, ";"));

        // print the rest of the parameters for the blocks on this read
        read_block_stats << "," << validBlocks.size() << "," << blocks.size() << std::endl;

        std::cout << std::endl;

        // Link blocks
        if (!validBlocks.empty()) {
            auto prev = validBlocks.cbegin();
            auto next = prev;
            ++next;
            for (; next != validBlocks.cend(); ++next) {
                if (std::abs(next->contigID) < std::abs(prev->contigID)) {
                    links.emplace_back(next->contigID * -1, prev->contigID * -1);
                    std::cout << links.back();
                } else if (std::abs(next->contigID) != std::abs(prev->contigID)) {
                    links.emplace_back(prev->contigID, next->contigID);
                    std::cout << links.back();
                }
                prev = next;
            }
        }
        std::cout << std::endl;
        read_link_stats << links.size() << std::endl;

        float pcReadCovered(matches.size() / (currentFileRecord.seq.size() - K + 1.0f));
        std::cout << "\nRead covered: " << matches.size() << "/" << currentFileRecord.seq.size() - K + 1 << std::endl;
        auto validBlockKmers(std::accumulate(validBlocks.begin(),validBlocks.end(), 0, [](uint32_t vbk, Block const &b){ return vbk+b.count; }));
        read_coverage_stats << matches.size() << "," << validBlockKmers << "," << currentFileRecord.seq.size()-K+1 << std::endl;
        std::cout << std::endl;

        for (const auto &l: links){
            read_contig_link << l.contig1 << ";" << l.contig2 << std::endl;
        }
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


    const std::unordered_set<KmerIDX> kmers;
    const std::unordered_map<int32_t, std::pair<uint64_t, uint32_t>> unitigs;

    char b2f[255];
    char b2r[255];

    // T = Total
    std::ofstream ctg_read;                     /// DEBUG: output only for debug mode
    mutable std::ofstream read_length;          /// INFO - Read, Length
    mutable std::ofstream read_coverage_stats;  /// INFO - Mapped kmers, ValidBlockKmers, TKmer
    mutable std::ofstream read_block_stats;     /// INFO - blockCtgDif, blockOffsetDif, Validblock, Tblock
    mutable std::ofstream read_link_stats;      /// INFO - # Links generated
};

#endif //SEQSORTER_CONTIGLINK_H_H
