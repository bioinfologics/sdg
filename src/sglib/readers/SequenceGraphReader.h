//
// Created by Luis Yanes (EI) on 17/11/2017.
//

#ifndef SG_SEQUENCEGRAPHREADER_H
#define SG_SEQUENCEGRAPHREADER_H
#include "Common.h"

/*
 * Pseudo-reader that gets its sequences from the nodes of the graph.
 */
class SequenceGraph;
struct GraphNodeReaderParams {
    uint32_t min_length;
    const SequenceGraph& sgp;
};

template<typename FileRecord>
class   GraphNodeReader {
public:
    explicit GraphNodeReader(GraphNodeReaderParams params) : params(params), numRecords(1), sg(params.sgp) {
        min_length=params.min_length;//TODO: use this
    }

    bool next_record(FileRecord& rec) {
        if (numRecords<sg.nodes.size()) {
            rec.id = numRecords;
            rec.seq = sg.nodes[numRecords].sequence;
            rec.name = std::to_string(numRecords);
            ++numRecords;
            stats.totalLength += rec.seq.size();
            return true;
        } else return false;
    }
    ReaderStats getSummaryStatistics() {
        stats.totalRecords = numRecords-1;
        return stats;
    }
private:
    const SequenceGraph& sg;
    GraphNodeReaderParams params;
    uint32_t numRecords;
    ReaderStats stats;
    uint32_t min_length;
};

#endif //SG_SEQUENCEGRAPHREADER_H
