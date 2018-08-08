//
// Created by Luis Yanes (EI) on 12/02/2018.
//

#ifndef BSG_LONGREADMAPPER_H
#define BSG_LONGREADMAPPER_H


#include <iostream>
#include <sglib/factories/KMerIDXFactory.h>
#include <sglib/datastores/LongReadsDatastore.hpp>
#include <sglib/SequenceGraph.hpp>

/**
 * TODO: Generate indices for LongReadMapping
 * READS FROM NODE -> IN: NODE ID OUT: READ set
 * NODES FROM READ -> IN: READ ID OUT: NODE set
 *
 * Stores the node, read, respective start and eds
 */
struct LongReadMapping {
    LongReadMapping() {}
    LongReadMapping(sgNodeID_t node, uint32_t read_id, int32_t nStart, int32_t nEnd, int32_t qStart, int32_t qEnd) :
    node(node), read_id(read_id), nStart(nStart), nEnd(nEnd), qStart(qStart), qEnd(qEnd) {}
    bool operator==(const LongReadMapping &other) const {
        return std::tie(node,read_id,nStart,nEnd,qStart,qEnd)
               == std::tie(other.node,other.read_id,other.nStart,other.nEnd,other.qStart,other.qEnd);
    }

    bool operator<(const LongReadMapping &other) const {
        return std::tie(node,read_id,nStart,nEnd,qStart,qEnd)
               < std::tie(other.node,other.read_id,other.nStart,other.nEnd,other.qStart,other.qEnd);
    }

    sgNodeID_t node = 0;        /// Node ID, sign represents direction
    uint32_t read_id = 0;       /// ID of the read from the Datastore   (this is never negative!)
    int32_t nStart = 0;         /// Position of the starting node kmer of this mapping
    int32_t nEnd = 0;           /// Position of the ending node kmer of this mapping
    int32_t qStart = 0;         /// Query start position
    int32_t qEnd = 0;           /// Query end position
    int32_t score = 0;          /// Alignment score
};

/**
 * Long read mapping to the graph, this class manages storage and computation of the alignments.
 */
class LongReadMapper {
    class StringKMerFactory : protected KMerFactory {
    public:
        explicit StringKMerFactory(uint8_t k) : KMerFactory(k) {}

        // TODO: Adjust for when K is larger than what fits in uint64_t!
        const bool create_kmers(const std::string &s, std::vector<std::pair<bool, uint64_t>> &mers) {
            fkmer=0;
            rkmer=0;
            last_unknown=0;
            uint64_t p(0);
            mers.reserve(s.size());
            while (p < s.size()) {
                //fkmer: grows from the right (LSB)
                //rkmer: grows from the left (MSB)
                fillKBuf(s[p], 0, fkmer, rkmer, last_unknown);
                p++;
                if (last_unknown >= K) {
                    if (fkmer <= rkmer) {
                        // Is fwd
                        mers.emplace_back(true, fkmer);
                    } else {
                        // Is bwd
                        mers.emplace_back(false, rkmer);
                    }
                }
            }
            return false;
        }
    };

    const SequenceGraph & sg;

    uint8_t k=15;
    int window_size = 750;
    int window_slide = window_size/3;
    int min_kmers = 11;
    int spread_pct =90;
    int max_num_score_nodes = 100;
    float min_spread = spread_pct/100.f;

    /**
     * Stores an index of the mappings of a node to all the mappings where it appears.
     * This index can be queried to get information about all reads that map to a node.
     */
    std::vector< std::vector < std::vector<LongReadMapping>::size_type > > reads_in_node;        /// Reads matching node

    void update_indexes_from_mappings();

public:

    LongReadMapper(SequenceGraph &sg, LongReadsDatastore &ds, uint8_t k=15);
    ~LongReadMapper();

    LongReadMapper operator=(const LongReadMapper &other);

    void update_graph_index();
    LongReadsDatastore& getLongReadsDatastore() {return datastore;}


    void map_reads(std::unordered_set<uint32_t> readIDs = {});

    void read(std::string filename);

    void read(std::ifstream &inf);

    void write(std::string filename);

    void write(std::ofstream &ofs);


    std::vector<KmerIDX> assembly_kmers;


    LongReadsDatastore datastore;
    /**
     * This public member stores a flat list of mappings from the reads, it is accessed using the mappings_in_node index
     * or the read_to_mappings index.
     */
    std::vector<LongReadMapping> mappings;
    /**
     * Stores an index of the resulting mappings of a single long read, for each long read, stores the position of it's mappings.
     * This index can be used to query all the nodes that map to a single read.
     */
    std::vector< std::vector < std::vector<LongReadMapping>::size_type > > read_to_mappings;    /// Nodes in the read, 0 or empty = unmapped
};


#endif //BSG_LONGREADMAPPER_H
