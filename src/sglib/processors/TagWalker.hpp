//
// Created by Bernardo Clavijo (EI) on 01/03/2018.
//

#ifndef BSG_TAGWALKER_HPP
#define BSG_TAGWALKER_HPP


#include <sglib/WorkSpace.hpp>

class TagWalker {
public:
    TagWalker(WorkSpace & _ws, std::pair<sgNodeID_t,sgNodeID_t> hspnp) : ws(_ws),nodeA(hspnp.first),nodeB(hspnp.second){
        tagsA=ws.linked_read_mappers[0].get_node_tags(nodeA);
        tagsB=ws.linked_read_mappers[0].get_node_tags(nodeB);
    };
    void remove_crosstalk();
    void dump_reads(std::string prefix);
    SequenceGraphPath walk(float min_winner,float max_looser);

private:
    WorkSpace &ws;
    sgNodeID_t nodeA;
    sgNodeID_t nodeB;
    std::unordered_set<bsg10xTag> tagsA;
    std::unordered_set<bsg10xTag> tagsB;
    std::unordered_set<bsg10xTag> tags_shared;
};


#endif //BSG_TAGWALKER_HPP
