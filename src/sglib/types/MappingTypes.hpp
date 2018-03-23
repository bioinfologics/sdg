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

typedef uint32_t prm10xTag_t;

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
        os << rm.node << "\t" << rm.unique_matches;
        return os;
    }

    sgNodeID_t node = 0;        /// Node ID
    uint64_t read_id = 0;       /// ID of the read from the Datastore
    int32_t first_pos = 0;      /// Position of the first node kmer of this mapping
    int32_t last_pos = 0;       /// Position of the last node kmer of this mapping
    int32_t unique_matches = 0; /// Number of unique kmer matches
    bool rev=false;
};

#endif //BSG_MAPPINGTYPES_HPP
