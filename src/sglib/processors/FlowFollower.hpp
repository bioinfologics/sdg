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
    void create_flows_all_fast();
    std::vector<std::unordered_set<uint64_t>> get_distinctive_kmers(std::vector<sgNodeID_t>);
    std::vector<std::unordered_set<uint64_t>> get_distinctive_kmers_truncated(std::vector<sgNodeID_t>);
    Flow flow_from_node(sgNodeID_t n,float min_winner,float max_looser);
    std::vector<SequenceGraphPath> skate_from_all(int min_node_flow, uint64_t min_path_length);
    SequenceGraphPath skate_from_node(sgNodeID_t);
    WorkSpace &ws;
    std::unordered_set<sgNodeID_t> nodes;
    std::map<sgNodeID_t,Flow> flows;

};


#endif //BSG_FLOWFOLLOWER_HPP
