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

const std::vector<LinkView> NodeView::next() const {
    auto fwl=dg->get_fw_links(id);
    std::vector<LinkView> r;
    r.reserve(fwl.size());
    for (auto &l:fwl) {
        r.emplace_back(dg->get_nodeview(l.dest),l.dist,l.support);
    }
    std::sort(r.begin(),r.end());
    return r;
}

const std::vector<LinkView> NodeView::prev() const {
    auto bwl=dg->get_bw_links(id);
    std::vector<LinkView> r;
    r.reserve(bwl.size());
    for (auto &l:bwl) {
        r.emplace_back(dg->get_nodeview(-l.dest),l.dist,l.support);
    }
    std::sort(r.begin(),r.end());
    return r;
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

std::vector<LinkedTag> NodeView::get_linked_tags(std::string datastore_name) const {
    std::vector<LinkedTag> tags;
    for (auto t:dg->sdg.ws.get_linked_reads_datastore(datastore_name).mapper.get_node_tags(llabs(id))) tags.push_back(t);
    return tags;
}

std::ostream &operator<<(std::ostream &os, const NodeView &nv) {
    os << "< NodeView: "<<nv.id<<" in "<<nv.dg->name<<" >";
    return os;
}

std::ostream &operator<<(std::ostream &os, const LinkView &ndv) {
    os << "< LinkView: "<<ndv.dist<<"bp to "<<ndv.node_view.node_id()<<" >";
    return os;
}