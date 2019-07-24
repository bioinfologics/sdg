//
// Created by Bernardo Clavijo (EI) on 2019-06-26.
//

#pragma once

#include <sdglib/graph/SequenceDistanceGraph.hpp>
#include <sdglib/types/GenericTypes.hpp>

class NodeDistanceView;

class NodeView {
public:
    NodeView(DistanceGraph * _dg,sgNodeID_t _n):dg(_dg),id(_n){};
    NodeView(NodeView const &o):dg(o.dg), id(o.id) {};
    friend std::ostream &operator<<(std::ostream &os, const NodeView &n);
    DistanceGraph graph() const;
    const std::string sequence() const;
    const uint64_t size() const;
    const std::vector<NodeDistanceView> next() const;
    const std::vector<NodeDistanceView> prev() const;
    const sgNodeID_t node_id() const {return sgNodeID_t(id);};
    std::vector<uint16_t> kmer_coverage(std::string kcovds_name,std::string kcovds_count_name) const;
    std::vector<uint16_t> kmer_coverage(int kcovds_idx,int kcovds_count_idx) const;
//    std::vector<seqID_t> preads(int prds_idx=0) const;
//    std::vector<seqID_t> preads(std::string prds_name) const;
//    std::vector<seqID_t> lireads(int lirds_idx=0) const;
//    std::vector<seqID_t> lireads(std::string lirds_name) const;
//    void li_tag_neihgbours(int lirds_idx=0) const;
//    void li_tag_neihgbours(std::string lirds_name) const;
//    void li_tags(int lirds_idx=0) const;
//    void li_tags(std::string lirds_name) const;
//    std::vector<seqID_t> loreads(int lords_idx=0) const;
//    std::vector<seqID_t> loreads(std::string lords_name) const;
//    std::vector<seqID_t> readpaths(int rpds_idx=0) const;
//    std::vector<seqID_t> readpaths(std::string rpds_name) const;

private:
    sgNodeID_t id;
    DistanceGraph * dg;
};

class NodeDistanceView {
public:
    NodeDistanceView(const NodeView &nv,const int32_t &d, const Support &s):node_view(nv),dist(d),sup(s){};
    NodeDistanceView(NodeDistanceView const &o):node_view(o.node_view),dist(o.dist),sup(o.sup) {};
    friend std::ostream &operator<<(std::ostream &os, const NodeDistanceView &ndv);
    const NodeView node() const {return NodeView(node_view);};
    const int32_t distance() const {return dist; };
    const Support support() const {return Support(sup);};
    bool operator<(const NodeDistanceView & other) const{
        if (dist<other.dist) return true;
        return node_view.node_id() < other.node_view.node_id();
    }
private:
    NodeView node_view;
    int32_t dist;
    Support sup;
};