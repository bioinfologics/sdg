//
// Created by Bernardo Clavijo (EI) on 2019-06-26.
//

#pragma once

#include <sdglib/graph/SequenceDistanceGraph.hpp>
#include <sdglib/types/GenericTypes.hpp>
#include <sdglib/mappers/LinkedReadsMapper.hpp>

class LinkView;

/**
 * @brief
 * The NodeView provides a read-only comforatable interface for graph traversal
 */
class NodeView {
public:
    NodeView(const DistanceGraph * const _dg,sgNodeID_t _n):dg(_dg),id(_n){};
    NodeView(NodeView const &o) = default;
    friend std::ostream &operator<<(std::ostream &os, const NodeView &n);
    std::string str() const{
        return "NodeView: Node "+std::to_string(id)+" in "+dg->name;
    }
    std::string repr() const{
        return "<NodeView: Node "+std::to_string(id)+" in "+dg->name+">";
    }
    const DistanceGraph & graph() const;
    bool operator==(const NodeView &o) const {return (id==o.id) and (dg==o.dg);};
    /**
     * @return Sequence of the underlying SequenceDistanceGraph node
     */
    std::string sequence() const;

    /**
     * @return Length of the node sequence in the graph
     */
    uint64_t size() const;

    /**
     * @return The reverse complement NodeView
     */
    NodeView rc() const;

    /**
     * @return A vector of LinkView for the next elements defined by the DistanceGraph
     */
    std::vector<LinkView> next() const;
    /**
     * @return A vector of LinkView for the previous elements defined by the DistanceGraph
     */
    std::vector<LinkView> prev() const;

    /**
     * @return A vector of NodeView for all nodes that have the same connections fw and bw than this node
     */
    std::vector<NodeView> parallels() const;

    /**
     * @return The id of the underlying SequenceDistanceGraph node
     */
    sgNodeID_t node_id() const {return sgNodeID_t(id);};

    /**
     * @brief Coverage of each kmer in the node
     * @param kcovds_name Name of the KmersCounter
     * @param kcovds_count_name Name of the count object in the KmersCounter
     * @return A vector with the kmer coverage count for each kmer in the node
     */
    std::vector<uint16_t> kmer_coverage(std::string kcovds_name,std::string kcovds_count_name) const;

    /**
     * @brief Coverage of each kmer in the node
     * @param kcovds_idx Index of the KmersCounter
     * @param kcovds_count_idx Index of the count object in the KmersCounter
     * @return A vector with the kmer coverage count for each kmer in this node
     */
    std::vector<uint16_t> kmer_coverage(int kcovds_idx,int kcovds_count_idx) const;

    float kci();
    /**
     * @brief Collect all paired reads from a PairedReadsDatastore referred to by name
     * @param datastore_name Name of the datastore to collect the reads from
     * @return Vector of all reads in the paired datastore that have mappings to this node
     */
    std::vector<seqID_t> get_paired_reads(std::string datastore_name) const;

    /**
     * @brief Collect all mappings of this node in a PairedReadsDatastore referred to by name
     * @param datastore_name Name of the datastore to collect the reads from
     * @return Vector of all mappings to this node in the PairedReadsMapper
     */
    std::vector<ReadMapping> get_paired_mappings(std::string datastore_name) const;
    /**
     * @brief Collect all Long reads from a LongReadsDatastore referred to by name
     * @param datastore_name Name of the datastore to collect the reads from
     * @return Vector of all reads in the paired datastore that have mappings to this node
     */
    std::vector<seqID_t> get_long_reads(std::string datastore_name) const;

    /**
     * @brief Collect all mappings of this node in a LongReadsDatastore referred to by name
     * @param datastore_name Name of the datastore to collect the reads from
     * @return Vector of all mappings to this node in the LongReadsMapper
     */
    std::vector<LongReadMapping> get_long_mappings(std::string datastore_name) const;
    /**
     * @brief Collect all Linked reads from a LinkedReadsDatastore referred to by name
     * @param datastore_name Name of the datastore to collect the reads from
     * @return Vector of all reads in the paired datastore that have mappings to this node
     */
    std::vector<seqID_t> get_linked_reads(std::string datastore_name) const;

    /**
     * @brief Collect all mappings of this node in a LinkedReadsDatastore referred to by name
     * @param datastore_name Name of the datastore to collect the reads from
     * @return Vector of all mappings to this node in the LinkedReadsMapper
     */
    std::vector<ReadMapping> get_linked_mappings(std::string datastore_name) const;

    /**
     * @brief Collect all tags of this node in a LinkedReadsDatastore referred to by name
     * @param datastore_name Name of the datastore to collect the reads from
     * @param Min number of reads mapped to the node to include a tag in the result
     * @return Vector of all tags that have reads on this node in the LinkedReadsMapper
     */
    std::vector<LinkedTag> get_linked_tags(std::string datastore_name, int min_read_count=3) const;

    /**
     * Vector of FW distances and neighbours from the nodeview, sorted in ascending distance (median if multi-link) order.
     * @param min_links
     * @return
     */
    std::vector<std::pair<int,sgNodeID_t>> fw_neighbours_by_distance( int min_links) const;

    std::vector<uint64_t> get_kmers(int K);

    std::unordered_set<uint64_t> get_linked_tags_kmers(std::string datastore_name, int K, int min_tag_cov);

    bool is_tip();
    bool is_bubble_side();
    bool is_canonical_repeat();

private:
    sgNodeID_t id;
    const DistanceGraph * dg;
};

/**
 * @brief
 * The LinkView provides a read-only interface to links in the DistanceGraph.
 * LinkView can only be generated using the next, prev methods in NodeView.
 *
 */
class LinkView {
public:
    LinkView(const NodeView &nv,const int32_t &d, const Support &s):node_view(nv),dist(d),sup(s){};
    LinkView(LinkView const &o) = default;
    friend std::ostream &operator<<(std::ostream &os, const LinkView &ndv);
    /**
     * @brief
     * NodeView of the DistanceGraph node referenced by the link
     * @return A NodeView of the traversed ID
     */
    const NodeView node() const {return NodeView(node_view);};

    /**
     * @brief
     * Distance to the neighbour node this LinkView was generated from
     * @return A negative value represents an overlap, a positive value represents a gap.
     */
    const int32_t distance() const {return dist; };

    /**
     * @brief Provides information about the origin of the Link
     * @return A summary of the information about the origin of the Link
     */
    const Support support() const {return Support(sup);};
    bool operator<(const LinkView & other) const{
        return std::make_tuple(dist, node_view.node_id()) < std::make_tuple(other.dist, other.node_view.node_id());
    }
private:
    NodeView node_view;
    int32_t dist;
    Support sup;
};