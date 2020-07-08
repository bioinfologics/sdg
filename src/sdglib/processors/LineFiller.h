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

    std::vector<sgNodeID_t> score_function(sgNodeID_t n1, sgNodeID_t n2, std::vector<SequenceDistanceGraphPath> paths);
    std::vector<sgNodeID_t> line_fill(std::vector<sgNodeID_t> anchor_path);

    WorkSpace & ws;
    LongReadsRecruiter& lrr;
};


#endif //SDG_LINEFILLER_H
