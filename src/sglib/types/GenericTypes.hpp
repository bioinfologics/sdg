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
typedef int32_t seqID_t; //first sequence is 0;

enum sgNodeStatus_t {sgNodeActive, sgNodeDeleted};

/**
 * The Node contains the sequence of a node and its status {Active, Deleted}
 */
class Node{
public:
    Node(std::string _seq, sgNodeStatus_t _status) : sequence(_seq), status(_status){};
    Node(std::string _seq) : sequence(_seq),status(sgNodeActive){};
    Node() = default;
    std::string sequence = "";
    sgNodeStatus_t status = sgNodeActive;
    bool is_canonical();
    void make_rc();

    friend std::ostream &operator<<(std::ostream &os, const Node &node) {
        if (node.sequence.length() > 20) {
            os << node.sequence.substr(0, 20) << " ... " << node.sequence.substr(node.sequence.length() - 20, 20);
        } else {
            os << node.sequence;
        }
        return os;
    }

};

/**
 * The Link represents a connection between one end of a sequence and the end of another sequence
 *
 * +AAAAAAAAA- connected to +BBBBBBBBBB-
 *
 * would be represented as the -A,B link. Links also contain a distance parameter in case of overlap the distance is
 * negative and in case of "scaffolding" link the distance would be positive.
 */
class Link{
public:
    Link(){};
    Link( sgNodeID_t _src, sgNodeID_t _dst, int32_t _dist) : source(_src), dest(_dst), dist(_dist) {};
    sgNodeID_t source = 0;
    sgNodeID_t dest = 0;
    int32_t dist = 0;

    bool operator==( const  Link);
    bool operator<(const Link)const;

    friend std::ostream &operator<<(std::ostream &os, const Link &link) {
        os << link.source << " -> " << link.dest;
        return os;
    }

};

/**
 * A node visitor contains the node ID, and distances in terms of NTs and Nodes from the starting node
 * This class is used as a helper for the depth_ and breath_fist_search functions
 */
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
