//
// Created by Bernardo Clavijo (EI) on 01/03/2018.
//

#ifndef BSG_TAGWALKER_HPP
#define BSG_TAGWALKER_HPP


#include <sglib/WorkSpace.hpp>

class TagWalker {
public:
    TagWalker(WorkSpace & _ws, std::pair<sgNodeID_t,sgNodeID_t> hspnp) :
            ws(_ws),nodeA(hspnp.first),nodeB(hspnp.second),pathA(_ws.sg,{}),pathB(_ws.sg,{}){
        tagsA=ws.linked_read_mappers[0].get_node_tags(nodeA);
        tagsB=ws.linked_read_mappers[0].get_node_tags(nodeB);
    };
    float remove_crosstalk();
    void dump_reads(std::string prefix);
    std::vector<SequenceGraphPath> walk(float min_winner,float max_looser);
    std::vector<std::unordered_set<uint64_t>> get_distinctive_kmers(std::vector<sgNodeID_t>);

private:
    WorkSpace &ws;
    sgNodeID_t nodeA;
    sgNodeID_t nodeB;
    std::unordered_set<bsg10xTag> tagsA;
    std::unordered_set<bsg10xTag> tagsB;
    std::unordered_set<bsg10xTag> tags_shared;
    SequenceGraphPath pathA,pathB;
};


#endif //BSG_TAGWALKER_HPP
