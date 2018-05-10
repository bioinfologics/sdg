//
// Created by Luis Yanes (EI) on 22/03/2018.
//

#ifndef BSG_GENERICTYPES_HPP
#define BSG_GENERICTYPES_HPP

#include <iostream>
#include <string>
#include <tuple>
#include <limits>
#include "hashing_helper.hpp"



typedef int64_t sgNodeID_t; //first node is 1; negatives are RC

enum sgNodeStatus_t {sgNodeActive,sgNodeDeleted};

class Node{
public:
    Node(std::string _seq, sgNodeStatus_t _status) : sequence(_seq), status(_status){};
    Node(std::string _seq) : sequence(_seq),status(sgNodeActive){};
    Node() = default;
    std::string sequence;
    sgNodeStatus_t status;
    bool is_canonical();
    void make_rc();
};

class Link{
public:
    Link( sgNodeID_t _src, sgNodeID_t _dst, int32_t _dist) : source(_src), dest(_dst), dist(_dist) {};
    sgNodeID_t source,dest;
    int32_t dist;

    bool operator==( const  Link);
    bool operator<(const Link)const;

    friend std::ostream &operator<<(std::ostream &os, const Link &link) {
        os << link.source << " -> " << link.dest;
        return os;
    }

};

struct nodeVisitor {
    sgNodeID_t node = 0;
    unsigned int dist = 0;
    unsigned int path_length = 0;
    nodeVisitor(sgNodeID_t n, unsigned int d, unsigned int p) : node(n), dist(d), path_length(p) {}
    nodeVisitor() = default;
    bool operator<(const nodeVisitor &o) const {return std::tie(node) < std::tie(o.node);}
    bool operator==(const nodeVisitor &o) const {return node == o.node;}
    nodeVisitor reverseDirection() const {
        return {-node, dist, path_length};
    }
    friend std::ostream &operator<<(std::ostream &os, const nodeVisitor &visitor) {
        os << visitor.node << ":" << visitor.path_length;
        return os;
    }
};

#endif //BSG_GENERICTYPES_HPP
