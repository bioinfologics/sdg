//
// Created by Bernardo Clavijo (EI) on 2019-06-26.
//

#include "NodeView.hpp"
#include <sdglib/graph/SequenceDistanceGraph.hpp>
#include <sdglib/workspace/WorkSpace.hpp>

DistanceGraph NodeView::graph() const {
    return *dg;
}

const std::string NodeView::sequence() const {
    return dg->sdg.get_node_sequence(id);
}

const uint64_t NodeView::size() const {
    return dg->sdg.get_node_size(id);
}

NodeView NodeView::rc() const {
    return dg->get_nodeview(-id);
}

const std::vector<LinkView> NodeView::next() const {
    auto fwl=dg->get_fw_links(id);
    std::vector<LinkView> r;
    r.reserve(fwl.size());
    for (auto &l:fwl) {
        r.emplace_back(dg->get_nodeview(l.dest),l.dist,l.support);
    }

    std::sort(r.begin(), r.end());
    return r;
}

const std::vector<LinkView> NodeView::prev() const {
    auto bwl=dg->get_bw_links(id);
    std::vector<LinkView> r;
    r.reserve(bwl.size());
    for (auto &l:bwl) {
        r.emplace_back(dg->get_nodeview(-l.dest),l.dist,l.support);
    }

    std::sort(r.begin(), r.end());
    return r;
}

const std::vector<NodeView> NodeView::parallels() const {
    auto p=prev(),n=next();
    if (p.empty() or n.empty()) return {};
    std::vector<NodeView> pars;
    std::set<sgNodeID_t> pnodes,nnodes;
    for (auto &plv:prev()) pnodes.insert(plv.node().node_id());
    for (auto &nlv:next()) nnodes.insert(nlv.node().node_id());
    for (auto &other:prev()[0].node().next()){
        if (other.node()==*this) continue;
        std::set<sgNodeID_t> opnodes,onnodes;
        for (auto &plv:other.node().prev()) opnodes.insert(plv.node().node_id());
        for (auto &nlv:other.node().next()) onnodes.insert(nlv.node().node_id());
        if (pnodes==opnodes and nnodes==onnodes) pars.push_back(other.node());
    }
    return pars;
}

std::vector<uint16_t> NodeView::kmer_coverage(std::string kcovds_name, std::string kcovds_count_name) const {
    return dg->sdg.ws.get_kmer_counter(kcovds_name).project_count(kcovds_count_name,sequence());
}

std::vector<uint16_t> NodeView::kmer_coverage(int kcovds_idx, int kcovds_count_idx) const {
    return dg->sdg.ws.kmer_counters[kcovds_count_idx].project_count(kcovds_count_idx,sequence());
}

std::vector<seqID_t> NodeView::get_paired_reads(std::string datastore_name) const {
    std::vector<seqID_t> reads;
    for(auto rm:dg->sdg.ws.get_paired_reads_datastore(datastore_name).mapper.reads_in_node[llabs(id)]) reads.push_back(rm.read_id);
    return reads;
}

std::vector<ReadMapping> NodeView::get_paired_mappings(std::string datastore_name) const {
    auto &mapper=dg->sdg.ws.get_paired_reads_datastore(datastore_name).mapper;
    auto mappings=dg->sdg.ws.get_paired_reads_datastore(datastore_name).mapper.reads_in_node[llabs(id)];
    if (id<0){
        for (auto &rm:mappings) {
            rm.node*=-1;
            rm.rev=not rm.rev;
            auto l=size()-rm.first_pos;
            auto f=size()-rm.last_pos;
            rm.first_pos=f;
            rm.last_pos=l;
        }
    }
    std::sort(mappings.begin(),mappings.end());
    return mappings;
}

std::vector<seqID_t> NodeView::get_long_reads(std::string datastore_name) const {
    std::vector<seqID_t> reads;
    for(auto rm:dg->sdg.ws.get_long_reads_datastore(datastore_name).mapper.reads_in_node[llabs(id)]) reads.push_back(rm);
    return reads;
}

std::vector<LongReadMapping> NodeView::get_long_mappings(std::string datastore_name) const {
    auto &mapper=dg->sdg.ws.get_long_reads_datastore(datastore_name).mapper;
    std::vector<LongReadMapping> mappings;
    for (auto lrid:get_long_reads(datastore_name)){
        for(auto &m:mapper.get_raw_mappings_from_read(lrid)){
            if (m.node==id) mappings.emplace_back(id,m.read_id,m.nStart,m.nEnd,m.qStart,m.qEnd,m.score);
            if (m.node==-id) mappings.emplace_back(id,m.read_id,size()-m.nStart,size()-m.nEnd,m.qStart,m.qEnd,m.score);
        }
    }
    std::sort(mappings.begin(),mappings.end());
    return mappings;
}

std::vector<seqID_t> NodeView::get_linked_reads(std::string datastore_name) const {
    std::vector<seqID_t> reads;
    for(auto rm:dg->sdg.ws.get_linked_reads_datastore(datastore_name).mapper.reads_in_node[llabs(id)]) reads.push_back(rm.read_id);
    return reads;
}

std::vector<ReadMapping> NodeView::get_linked_mappings(std::string datastore_name) const {
    auto &mapper=dg->sdg.ws.get_linked_reads_datastore(datastore_name).mapper;
    auto mappings=dg->sdg.ws.get_linked_reads_datastore(datastore_name).mapper.reads_in_node[llabs(id)];
    if (id<0){
        for (auto &rm:mappings) {
            rm.node*=-1;
            rm.rev=not rm.rev;
            auto l=size()-rm.first_pos;
            auto f=size()-rm.last_pos;
            rm.first_pos=f;
            rm.last_pos=l;
        }
    }
    std::sort(mappings.begin(),mappings.end());
    return mappings;
}

std::vector<LinkedTag> NodeView::get_linked_tags(std::string datastore_name, int min_read_count) const {

    std::map<LinkedTag, uint32_t> tag_counts;
    for (auto mapping: dg->sdg.ws.get_linked_reads_datastore(datastore_name).mapper.reads_in_node[llabs(id)]){
        auto t = dg->sdg.ws.get_linked_reads_datastore(datastore_name).get_read_tag(mapping.read_id);
        tag_counts[t]++;
    }

    std::vector<LinkedTag> tags;
    for (auto &t: tag_counts) {
        if (t.second > min_read_count) {
            tags.push_back(t.first);
        }
    }
    return tags;
}

std::ostream &operator<<(std::ostream &os, const NodeView &nv) {
    os << "NodeView: Node "<<nv.id<<" in "<<nv.dg->name;
    return os;
}

std::ostream &operator<<(std::ostream &os, const LinkView &ndv) {
    os << "LinkView: "<<ndv.dist<<"bp to Node "<<ndv.node_view.node_id();
    return os;
}

std::vector<std::pair<int,sgNodeID_t>> NodeView::fw_neighbours_by_distance(int min_links) const {
    return dg->fw_neighbours_by_distance(id, min_links);
}

std::vector<uint64_t> NodeView::get_kmers(int K){
    std::vector<uint64_t> node_kmers;
    auto kf = StringKMerFactory(K);
    kf.create_kmers(sequence(), node_kmers);
    return node_kmers;
};

std::unordered_set<uint64_t> NodeView::get_linked_tags_kmers(std::string datastore_name, int K, int min_tag_cov){
    // Get tags
//    auto tags = get_linked_reads(datastore_name);
    std::set<LinkedTag> tags = dg->sdg.ws.get_linked_reads_datastore(datastore_name).mapper.get_node_tags(llabs(id));

    // Get kmers from tags
    ReadSequenceBuffer ds_buffer(dg->sdg.ws.get_linked_reads_datastore(datastore_name));
    return dg->sdg.ws.get_linked_reads_datastore(datastore_name).get_tags_kmers(K, min_tag_cov, tags, ds_buffer);
};