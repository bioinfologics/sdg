//
// Created by Bernardo Clavijo (EI) on 2019-06-26.
//

#pragma once

#include <sdglib/graph/SequenceDistanceGraph.hpp>

class NodeDistanceView;

class NodeView {
public:
    NodeView(DistanceGraph &_dg,sgNodeID_t _n):dg(_dg),node_id(_n){};
    NodeView& operator=(NodeView const &o);
    friend std::ostream &operator<<(std::ostream &os, const NodeView &n);
    std::string sequence() const;
    uint64_t size() const;
    std::vector<NodeDistanceView> next() const;
    std::vector<NodeDistanceView> prev() const;
//    std::vector<uint16_t> kmer_coverage(std::string kcovds_name,int kcovds_count_name) const;
//    std::vector<uint16_t> kmer_coverage(int kcovds_idx,int kcovds_count_idx) const;
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
    DistanceGraph &dg;
    sgNodeID_t node_id;
};

class NodeDistanceView {
public:
    NodeDistanceView(const NodeView &nv,const int32_t &d, const Support &s):node(nv),distance(d),support(s){};
    NodeDistanceView& operator=(NodeDistanceView const &o);
    friend std::ostream &operator<<(std::ostream &os, const NodeDistanceView &ndv);
    NodeView node;
    int32_t distance;
    Support support;
};