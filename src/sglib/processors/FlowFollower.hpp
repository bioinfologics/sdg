//
// Created by Bernardo Clavijo (EI) on 08/03/2018.
//

#ifndef BSG_FLOWFOLLOWER_HPP
#define BSG_FLOWFOLLOWER_HPP

#include <sglib/WorkSpace.hpp>

class Flow {
public:
    std::vector<sgNodeID_t> nodes;
};

/**
 * @brief The FlowFollower creates Flows from every node, and then finds paths that follow the flows through the graph
 */
class FlowFollower {
public:
    FlowFollower(WorkSpace & _ws,std::unordered_set<sgNodeID_t> _nodes={}):ws(_ws),nodes(_nodes){};
    void create_flows();
    std::vector<std::unordered_set<uint64_t>> get_distinctive_kmers(std::vector<sgNodeID_t>);
    std::vector<std::unordered_set<uint64_t>> get_distinctive_kmers_truncated(std::vector<sgNodeID_t>);
    Flow flow_from_node(sgNodeID_t n,float min_winner,float max_looser);
    WorkSpace &ws;
    std::unordered_set<sgNodeID_t> nodes;
    std::map<sgNodeID_t,Flow> flows;

};


#endif //BSG_FLOWFOLLOWER_HPP
