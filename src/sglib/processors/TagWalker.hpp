//
// Created by Bernardo Clavijo (EI) on 01/03/2018.
//

#ifndef BSG_TAGWALKER_HPP
#define BSG_TAGWALKER_HPP


#include <sglib/workspace/WorkSpace.hpp>

class TagWalker {
public:
    TagWalker(WorkSpace & _ws, std::pair<sgNodeID_t,sgNodeID_t> hspnp) :
            ws(_ws),nodeA(hspnp.first),nodeB(hspnp.second),pathA(_ws.getGraph(),{}),pathB(_ws.getGraph(),{}){
        tagsA=ws.getLinkedReadMappers()[0].get_node_tags(nodeA);
        tagsB=ws.getLinkedReadMappers()[0].get_node_tags(nodeB);
    };
    float remove_crosstalk();
    void dump_reads(std::string prefix);
    std::vector<SequenceGraphPath> walk(float min_winner,float max_looser);
    std::vector<std::unordered_set<uint64_t>> get_distinctive_kmers(std::vector<sgNodeID_t>);

private:
    WorkSpace &ws;
    sgNodeID_t nodeA;
    sgNodeID_t nodeB;
    std::set<bsg10xTag> tagsA;
    std::set<bsg10xTag> tagsB;
    std::set<bsg10xTag> tags_shared;
    SequenceGraphPath pathA,pathB;
};


#endif //BSG_TAGWALKER_HPP
