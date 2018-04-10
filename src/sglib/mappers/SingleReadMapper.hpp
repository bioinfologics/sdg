//
// Created by Luis Yanes (EI) on 09/04/2018.
//

#ifndef BSG_SINGLEREADMAPPER_HPP
#define BSG_SINGLEREADMAPPER_HPP


#include <sglib/types/GenericTypes.hpp>
#include <cstdint>
#include <vector>
#include <sglib/types/MappingTypes.hpp>
#include <sglib/types/KmerTypes.hpp>
#include <sglib/graph/SequenceGraph.hpp>
#include <sglib/datastores/LongReadsDatastore.hpp>
#include <sglib/indexers/minSketchIndex.hpp>

class SingleReadMapper {
    SequenceGraph & sg;
    LongReadsDatastore &datastore;
    std::unordered_map<uint64_t, graphPosition> kmer_to_graphposition;
    uint64_t memlimit;
    uint k;
    std::vector<std::vector<ReadMapping>> reads_in_node;
    std::vector<sgNodeID_t> read_to_node;//id of the main node if mapped, set to 0 to remap on next process
    std::vector<LongReadMapping> mappings;
    minSketchIndex index;

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
    SingleReadMapper(uint k, SequenceGraph &_sg, LongReadsDatastore &_datastore) : k(k), sg(_sg),datastore(_datastore), index(_sg, k, 1){
        reads_in_node.resize(sg.nodes.size());
    };

    void map_reads(std::unordered_set<uint64_t> const &  reads_to_remap={});

    std::vector<SingleReadMapper::MatchOffset> getMatchOffsets(std::string &query) {
        std::vector<MatchOffset> matches;
        uint32_t sketch_not_in_index(0);
        uint32_t sketch_in_index(0);
        StrandedMinimiserSketchFactory kf(k, 1);
        std::unordered_set<MinPosIDX> sketch;
        sketch.reserve(query.size());

        auto read_sketch_num(kf.getMinSketch(query, sketch));

        for (auto sk : sketch) {
            auto foundKey(index.find(sk.hash));
            if (foundKey == index.end()) {
                sketch_not_in_index++;
                continue;
            }
            sketch_in_index++;
            if (foundKey->second.size() < index.get_avg_nodes_per_kmer()) {
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

    std::vector<ReadMapping>
    map_read(uint32_t readID, std::string sequence, std::ofstream &matchOutput, std::ofstream &blockOutput) {
        std::vector<ReadMapping> paths;

        std::vector<unsigned int> nodeMatchCount(sg.nodes.size());

        auto matches = getMatchOffsets(sequence);
        auto sortedByContigRead = matches;

        std::for_each(matches.begin(),matches.end(),
                      [&nodeMatchCount](const MatchOffset &m) { nodeMatchCount[std::abs(m.dirContig)]+=1; } );

        std::sort(matches.begin(),matches.end(), MatchOffset::byReadPos());
        std::sort(sortedByContigRead.begin(), sortedByContigRead.end(), MatchOffset::byReadContig());

        for (const auto &m:matches) {
            matchOutput << "@MATCH " << readID << ","
                        << sg.nodes[std::abs(m.dirContig)].sequence.length() << "," << m
                        << std::endl;
        }

        return paths;
    }

};


#endif //BSG_SINGLEREADMAPPER_HPP
