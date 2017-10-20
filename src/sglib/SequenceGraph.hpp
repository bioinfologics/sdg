//
// Created by Bernardo Clavijo (EI) on 18/10/2017.
//

#ifndef SG_SEQUENCEGRAPH_HPP
#define SG_SEQUENCEGRAPH_HPP

#include <vector>
#include <string>
#include <unordered_map>

typedef int64_t sgNodeID_t; //first node is 1; negatives are RC

class SequenceGraph;
class Node{
public:
    Node(std::string _seq) : sequence(_seq){};
    std::string sequence;
    uint8_t status;
    bool is_canonical();
    void make_rc();
    //std::shared_ptr<SequenceGraph> sub_sg; <--- someday
};

class Link{
public:
    Link( sgNodeID_t _src, sgNodeID_t _dst, int32_t _dist) : source(_src), dest(_dst), dist(_dist) {};
    sgNodeID_t source,dest;
    int32_t dist;
};

class SequenceGraph {
public:
    SequenceGraph(){};
    void load_from_gfa(std::string filename);
    void write_to_gfa(std::string filename);
    sgNodeID_t add_node(Node n);

private:
    std::vector<Node> nodes;
    std::vector<std::vector<Link>> links;

};

class SequenceGraphPath {
public:
    std::vector<sgNodeID_t> nodes;
    explicit SequenceGraphPath(std::shared_ptr<SequenceGraph> _sg, std::vector<sgNodeID_t> _nodes={})  : sg(_sg) ,nodes(_nodes) {};
private:
    std::shared_ptr<SequenceGraph> sg;
};



#endif //SG_SEQUENCEGRAPH_HPP
