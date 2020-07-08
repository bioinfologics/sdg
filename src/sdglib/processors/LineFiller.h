//
// Created by Gonzalo Garcia (EI) on 2020-07-07.
//

#ifndef SDG_LINEFILLER_H
#define SDG_LINEFILLER_H


#include <vector>
#include <sdglib/graph/SequenceDistanceGraph.hpp>
#include <sdglib/mappers/LongReadsRecruiter.hpp>
#include <sdglib/views/NodeView.hpp>

class LineFiller {
public:
    LineFiller(WorkSpace & _ws, LongReadsRecruiter& _lrr):ws(_ws), lrr(_lrr){};

    uint32_t score_function(std::vector<sgNodeID_t> path);
    std::vector<sgNodeID_t> line_fill(std::vector<sgNodeID_t> anchor_path);


    void populate_matches(const LongReadsRecruiter& lrr);

    WorkSpace & ws;
    LongReadsRecruiter& lrr;
    std::vector<int> node_matches;
};


#endif //SDG_LINEFILLER_H
