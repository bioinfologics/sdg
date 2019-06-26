//
// Created by Bernardo Clavijo (EI) on 2019-06-26.
//

#pragma once

#include <sdglib/graph/DistanceGraph.hpp>

class NodeDistanceView;

class NodeView {
    NodeView(const DistanceGraph &_dg,sgNodeID_t _n):dg(_dg),node_id(_n){};
    std::vector<NodeDistanceView> next() const;
    std::vector<NodeDistanceView> prev() const;
    std::string sequence() const;
    std::vector<uint16_t> kmer_coverage(std::string kcovds_name,int kcovds_count_name) const;
    std::vector<uint16_t> kmer_coverage(int kcovds_idx,int kcovds_count_idx) const;
    std::vector<seqID_t> preads(int prds_idx=0) const;
    std::vector<seqID_t> preads(std::string prds_name) const;
    std::vector<seqID_t> lireads(int lirds_idx=0) const;
    std::vector<seqID_t> lireads(std::string lirds_name) const;
    void li_tag_neihgbours(int lirds_idx=0) const;
    void li_tag_neihgbours(std::string lirds_name) const;
    void li_tags(int lirds_idx=0) const;
    void li_tags(std::string lirds_name) const;
    std::vector<seqID_t> loreads(int lords_idx=0) const;
    std::vector<seqID_t> loreads(std::string lords_name) const;
    std::vector<seqID_t> readpaths(int rpds_idx=0) const;
    std::vector<seqID_t> readpaths(std::string rpds_name) const;
    const DistanceGraph &dg;
    const sgNodeID_t node_id;
};

class NodeDistanceView {
    NodeView node;
    int32_t distance;
    Support support;
};