//
// Created by Bernardo Clavijo (EI) on 2019-06-26.
//

#pragma once

#include <sdglib/graph/SequenceDistanceGraph.hpp>
#include <sdglib/types/GenericTypes.hpp>
#include <sdglib/mappers/LinkedReadsMapper.hpp>

class LinkView;

class NodeView {
public:
    NodeView(DistanceGraph * _dg,sgNodeID_t _n):dg(_dg),id(_n){};
    NodeView(NodeView const &o):dg(o.dg), id(o.id) {};
    friend std::ostream &operator<<(std::ostream &os, const NodeView &n);
    DistanceGraph graph() const;
    const std::string sequence() const;
    const uint64_t size() const;
    const std::vector<LinkView> next() const;
    const std::vector<LinkView> prev() const;
    const sgNodeID_t node_id() const {return sgNodeID_t(id);};
    std::vector<uint16_t> kmer_coverage(std::string kcovds_name,std::string kcovds_count_name) const;
    std::vector<uint16_t> kmer_coverage(int kcovds_idx,int kcovds_count_idx) const;
    std::vector<seqID_t> get_paired_reads(std::string datastore_name) const;
    std::vector<ReadMapping> get_paired_mappings(std::string datastore_name) const;
    std::vector<seqID_t> get_long_reads(std::string datastore_name) const;
    std::vector<LongReadMapping> get_long_mappings(std::string datastore_name) const;
    std::vector<seqID_t> get_linked_reads(std::string datastore_name) const;
    std::vector<ReadMapping> get_linked_mappings(std::string datastore_name) const;
    std::vector<LinkedTag> get_linked_tags(std::string datastore_name) const;

private:
    sgNodeID_t id;
    DistanceGraph * dg;
};

class LinkView {
public:
    LinkView(const NodeView &nv,const int32_t &d, const Support &s):node_view(nv),dist(d),sup(s){};
    LinkView(LinkView const &o):node_view(o.node_view),dist(o.dist),sup(o.sup) {};
    friend std::ostream &operator<<(std::ostream &os, const LinkView &ndv);
    const NodeView node() const {return NodeView(node_view);};
    const int32_t distance() const {return dist; };
    const Support support() const {return Support(sup);};
    bool operator<(const LinkView & other) const{
        if (dist<other.dist) return true;
        return node_view.node_id() < other.node_view.node_id();
    }
private:
    NodeView node_view;
    int32_t dist;
    Support sup;
};