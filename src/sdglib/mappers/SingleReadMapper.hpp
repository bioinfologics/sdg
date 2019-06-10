//
// Created by Luis Yanes (EI) on 09/04/2018.
//

#ifndef BSG_SINGLEREADMAPPER_HPP
#define BSG_SINGLEREADMAPPER_HPP


#include <sdglib/types/GenericTypes.hpp>
#include <cstdint>
#include <vector>
#include <sdglib/types/MappingTypes.hpp>
#include <sdglib/types/KmerTypes.hpp>
#include <sdglib/graph/SequenceDistanceGraph.hpp>
#include <sdglib/datastores/LongReadsDatastore.hpp>

class SingleReadMapper {
    SequenceGraph & sg;
    LongReadsDatastore &datastore;
    std::unordered_map<uint64_t, graphPosition> kmer_to_graphposition;
    uint64_t memlimit;
    unsigned int k;
    std::vector<std::vector<ReadMapping>> reads_in_node;
    std::vector<sgNodeID_t> read_to_node;//id of the main node if mapped, set to 0 to remap on next process
    std::vector<LongReadMapping> mappings;

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
                auto A = std::abs(a.dirContig);
                auto B = std::abs(b.dirContig);
                return std::tie(A, a.readPos) < std::tie(B, b.readPos);
            }
        };

        friend class byReadContig;
        struct byReadContig {
            bool operator()(const MatchOffset &a, const MatchOffset &b) {
                auto A = std::abs(a.dirContig);
                auto B = std::abs(b.dirContig);
                return std::tie(a.readPos,A) < std::tie(b.readPos, B);
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

public:
    SingleReadMapper(unsigned int k, SequenceGraph &_sg, LongReadsDatastore &_datastore) : k(k), sg(_sg),datastore(_datastore) {
        reads_in_node.resize(sg.nodes.size());
    };

    void map_reads(std::unordered_set<uint64_t> const &  reads_to_remap={});

};


#endif //BSG_SINGLEREADMAPPER_HPP
