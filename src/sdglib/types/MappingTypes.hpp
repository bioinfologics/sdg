//
// Created by Luis Yanes (EI) on 23/03/2018.
//

#ifndef BSG_MAPPINGTYPES_HPP
#define BSG_MAPPINGTYPES_HPP

#include <cstdint>
#include <strings.h>
#include <ostream>
#include <tuple>
#include "GenericTypes.hpp"
#include "hashing_helper.hpp"

class ReadMapping {
public:
    ReadMapping(){
        //just clean the structure, OSX doesn't give you clean memory
        bzero(this, sizeof(ReadMapping));
    }
    bool operator==(const ReadMapping &other){
        return this==&other;
    };
    bool operator<(const ReadMapping &other) const {
        return std::tie(node, read_id) < std::tie(other.node, other.read_id);
    };
    void merge(const ReadMapping &other){};
    friend std::ostream& operator<<(std::ostream& os, const ReadMapping& rm) {
        os << "ReadMapping " << rm.node << "\t" << rm.unique_matches;
        return os;
    }

    sgNodeID_t node = 0;        /// Node ID
    uint64_t read_id = 0;       /// ID of the read from the Datastore
    int32_t first_pos = 0;      /// Position of the first node kmer of this mapping
    int32_t last_pos = 0;       /// Position of the last node kmer of this mapping
    int32_t unique_matches = 0; /// Number of unique kmer matches
    bool rev=false;


};
struct ReadMapping_hash {
    size_t operator()(const ReadMapping& lr) const {
        std::tuple<sgNodeID_t , uint64_t , int32_t > tp (lr.node,lr.read_id,lr.first_pos);
        sdglib::hash<std::tuple<sgNodeID_t , uint64_t , int32_t>> h;
        return h (tp);
    }
};

/**
 * READS FROM NODE -> IN: NODE ID OUT: READ set
 * NODES FROM READ -> IN: READ ID OUT: NODE set
 *
 * Stores the node, read, respective start and eds
 */
struct LongReadMapping {
    LongReadMapping() {}
    LongReadMapping(sgNodeID_t nodeID, uint32_t readID, int32_t nstart = 0, int32_t nend = 0, int32_t qstart = 0, int32_t qend = 0, int32_t score = 0) :
    node(nodeID), read_id(readID), nStart(nstart), nEnd(nend), qStart(qstart), qEnd(qend), score(score) {}

    bool operator==(const LongReadMapping &other) const {
        return std::tie(node,read_id,nStart,nEnd,qStart,qEnd)
               == std::tie(other.node,other.read_id,other.nStart,other.nEnd,other.qStart,other.qEnd);
    }

    bool operator!=(const LongReadMapping &other) const {
        return !(*this==other);
    }

    bool operator<(const LongReadMapping &other) const {
        return std::tie(read_id,qStart,qEnd,node,nStart,nEnd)
               < std::tie(other.read_id,other.qStart,other.qEnd,other.node,other.nStart,other.nEnd);
    }

    friend std::ostream& operator<<(std::ostream& os, const LongReadMapping& m){
        os << "LongReadMapping " << m.read_id << " (" << m.qStart << ":" << m.qEnd << ") -> "
                                  << m.node << " (" << m.nStart << ":" << m.nEnd << ")  " <<m.score<< " hits";
        return os;
    }

    sgNodeID_t node = 0;        ///< Node ID, sign represents direction
    uint32_t read_id = 0;       ///< ID of the read from the Datastore (this is never negative!)
    int32_t nStart = 0;         ///< Position of the starting node kmer of this mapping
    int32_t nEnd = 0;           ///< Position of the ending node kmer of this mapping
    int32_t qStart = 0;         ///< Query start position
    int32_t qEnd = 0;           ///< Query end position
    int32_t score = 0;          ///< Alignment score

};

struct LongReadMapping_hash {
    size_t operator()(const LongReadMapping& lr) const {
        std::tuple<sgNodeID_t , uint64_t , int32_t > tp (lr.node,lr.read_id,lr.nStart);
        sdglib::hash<std::tuple<sgNodeID_t , uint64_t , int32_t>> h;
        return h (tp);
    }
};


struct SequenceMatch {
    SequenceMatch() {}
    SequenceMatch(sgNodeID_t seq_id1,sgNodeID_t seq_id2,int32_t start1,int32_t end1,int32_t start2,int32_t end2,int32_t hits) :
            seq_id1(seq_id1), seq_id2(seq_id2), start1(start1), end1(end1), start2(start2), end2(end2), hits(hits) {}

    bool operator==(const SequenceMatch &other) const {
        return std::tie(seq_id1,seq_id2,start1,end1,start2,end2,hits)
               == std::tie(other.seq_id1,other.seq_id2,other.start1,other.end1,other.start2,other.end2,other.hits);
    }

    bool operator!=(const SequenceMatch &other) const {
        return !(*this==other);
    }

    bool operator<(const SequenceMatch &other) const {
        return std::tie(seq_id1,seq_id2,start1,end1,start2,end2,hits)
               < std::tie(other.seq_id1,other.seq_id2,other.start1,other.end1,other.start2,other.end2,other.hits);
    }

    friend std::ostream& operator<<(std::ostream& os, const SequenceMatch& m){
        os << "SequenceMatch " << m.seq_id1 << " (" << m.start1 << ":" << m.end1 << ") <-> "
           << m.seq_id2 << " (" << m.start2 << ":" << m.end2 << ")  " <<m.hits<< " hits";
        return os;
    }

    sgNodeID_t seq_id1 = 0;        ///< Node ID, sign represents direction
    sgNodeID_t seq_id2 = 0;       ///< ID of the read from the Datastore (this is never negative!)
    int32_t start1 = 0;         ///< Position of the starting node kmer of this mapping
    int32_t end1 = 0;           ///< Position of the ending node kmer of this mapping
    int32_t start2 = 0;         ///< Query start position
    int32_t end2 = 0;           ///< Query end position
    int32_t hits = 0;          ///< Alignment score

};

struct SequenceMatch_hash {
    size_t operator()(const SequenceMatch& lr) const {
        std::tuple<sgNodeID_t , uint64_t , int32_t > tp (lr.seq_id1,lr.seq_id2,lr.start1);
        sdglib::hash<std::tuple<sgNodeID_t , uint64_t , int32_t>> h;
        return h (tp);
    }
};

#endif //BSG_MAPPINGTYPES_HPP
