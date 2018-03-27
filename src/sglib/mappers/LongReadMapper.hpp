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



class LongReadMapper {
    template <typename A, typename B>
    std::multimap<B, A> flip_map(std::map<A,B> & src) {

        std::multimap<B,A> dst;

        for(typename std::map<A, B>::const_iterator it = src.begin(); it != src.end(); ++it)
            dst.insert(std::pair<B, A>(it -> second, it -> first));

        return dst;
    }

    SequenceGraph & sg;
    minSketchIndex index;
    uint8_t k=15;
    uint8_t w=5;

    struct MatchOffset{
        MatchOffset(){};
        MatchOffset(sgNodeID_t dirContig, uint32_t readPos, uint32_t offset) :
                dirContig(dirContig), readPos(readPos), offset(offset) {};

        friend std::ostream &operator<<(std::ostream &os, const MatchOffset &match) {
            os << "" << match.dirContig << ","
               << match.readPos << ","
               << match.offset
               << "";
            return os;
        }
        sgNodeID_t dirContig;  // Sign indicates direction
        uint32_t readPos;   // Pos on the read
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

    std::vector<Match> getMatches(std::string &seq);
    std::vector<MatchOffset> getMatchOffsets(std::string &query);

public:
    LongReadMapper(uint8_t k, uint8_t w, SequenceGraph &sg) : sg(sg), k(k), w(w), index(sg, k, w) {}
    std::vector<LongReadMapping> map_read(FastqRecord read, std::ofstream &matchOutput, std::ofstream &blockOutput);
    uint64_t map_reads2(std::string &filename, uint32_t error);
    int32_t getWinner(std::multimap<uint32_t , int32_t> ranking, uint min_window_matches, float min_match_spread);
    std::vector<LongReadMapping> map_reads(std::string &filename);
};


#endif //BSG_LONGREADMAPPER_H
