//
// Created by Ben Ward (EI) on 02/10/2018.
//

#include "jlcxx/jlcxx.hpp"
#include "jlcxx/const_array.hpp"
#include <sdglib/graph/SequenceDistanceGraph.hpp>
#include <sdglib/workspace/WorkSpace.hpp>
#include <sdglib/processors/LinkageUntangler.hpp>
#include <sdglib/mappers/PairedReadsMapper.hpp>
#include <sdglib/graph/DistanceGraph.hpp>

namespace jlcxx {
    template<> struct IsBits<sgNodeStatus_t> : std::true_type {};
}

JLCXX_MODULE define_julia_module(jlcxx::Module& mod)
{
    // Basic Types //
/*
    // Read mapping types
    mod.add_type<ReadMapping>("ReadMapping")
            .constructor()
            .method("merge", &ReadMapping::merge);

    mod.add_type<LongReadMapping>("LongReadMapping")
            .constructor();
*/
    // Kmer types
    mod.add_type<KmerIDX>("KmerIDX")
            .constructor()
            .constructor<uint64_t>()
            .constructor<uint64_t, int32_t, uint32_t, uint8_t>()
            .method("isless", &KmerIDX::operator<)
            .method("isgreater", &KmerIDX::operator>)
            .method("isequal", &KmerIDX::operator==)
            .method("merge", &KmerIDX::merge)
            .method("kmer", [](KmerIDX& kidx) -> const uint64_t { return kidx.kmer; })
            .method("contigID", [](KmerIDX& kidx) -> const int32_t { return kidx.contigID; })
            .method("pos", [](KmerIDX& kidx) -> const uint32_t { return kidx.pos; })
            .method("count", [](KmerIDX& kidx) -> const uint8_t { return kidx.count; })
            .method("<<", [](KmerIDX& kidx){ std::cout << kidx; });

    mod.add_type<KmerIDX128>("KmerIDX128")
            .constructor()
            //.constructor<__uint128_t>()
            //.constructor<__uint128_t, int32_t, uint32_t, uint8_t>()
            .method("isless", &KmerIDX128::operator<)
            .method("isgreater", &KmerIDX128::operator>)
            .method("isequal", &KmerIDX128::operator==)
            .method("merge", &KmerIDX128::merge)
            .method("kmer", [](KmerIDX128& kidx) -> const __uint128_t { return kidx.kmer; })
            .method("contigID", [](KmerIDX128& kidx) -> const int32_t { return kidx.contigID; })
            .method("pos", [](KmerIDX128& kidx) -> const uint32_t { return kidx.pos; })
            .method("count", [](KmerIDX128& kidx) -> const uint8_t { return kidx.count; });

    mod.add_type<graphStrandPos>("graphStrandPos")
            .constructor()
            .constructor<sgNodeID_t, uint32_t>()
            .method("isequal", &graphStrandPos::operator==)
            .method("node", [](graphStrandPos& gsp){ return gsp.node; })
            .method("pos", [](graphStrandPos& gsp){ return gsp.pos; });

    mod.add_type<MinPosIDX>("MinPosIDX")
            .constructor()
            .constructor<uint64_t, int32_t>()
            .method("isless", &MinPosIDX::operator<)
            .method("isequal", &MinPosIDX::operator==)
            .method("hash", [](MinPosIDX& mp){ return mp.hash; })
            .method("pos", [](MinPosIDX& mp){ return mp.pos; })
            .method("<<", [](MinPosIDX& mp) { std::cout << mp; });

    mod.add_type<kmerPos>("kmerPos")
            .constructor()
            .constructor<uint64_t, uint32_t, int32_t>()
            .method("isless", &kmerPos::operator<)
            .method("kmer", [](kmerPos& kp) -> const uint64_t { return kp.kmer; })
            .method("contigID", [](kmerPos& kp) -> const int32_t { return kp.contigID; })
            .method("offset", [](kmerPos& kp) -> const int32_t { return kp.offset; });

    mod.add_type<nodeVisitor>("nodeVisitor")
            .constructor()
            .constructor<sgNodeID_t, unsigned int, unsigned int>()
            .method("reverseDirection", &nodeVisitor::reverseDirection)
            .method("<<", [](nodeVisitor& nv){ std::cout << nv; })
            .method("isless", &nodeVisitor::operator<)
            .method("isequal", &nodeVisitor::operator==)
            .method("node", [](nodeVisitor& nv) -> const sgNodeID_t { return nv.node; })
            .method("dist", [](nodeVisitor& nv) -> const unsigned int { return nv.dist; })
            .method("path_length", [](nodeVisitor& nv) -> const unsigned int { return nv.path_length; });

    mod.add_bits<sgNodeStatus_t>("sgNodeStatus_t");
    mod.set_const("sgNodeActive", sgNodeActive);
    mod.set_const("sgNodeDeleted", sgNodeDeleted);

    mod.add_type<Node>("Node")
            .constructor<std::string, sgNodeStatus_t>()
            .constructor<std::string>()
            .method("is_canonical", &Node::is_canonical)
            .method("make_rc", &Node::make_rc)
            .method("sequence", [](Node& nd){ return nd.sequence; })
            .method("status", [](Node& nd){ return nd.status; })
            .method("==", &Node::operator==)
            .method("<<", [](Node& n){ std::cout << n; });

    mod.add_type<Link>("Link")
            .constructor()
            .constructor<sgNodeID_t, sgNodeID_t, int32_t>()
            .method("source", [](Link& lnk){ return lnk.source; })
            .method("dest", [](Link& lnk){ return lnk.dest; })
            .method("dist", [](Link& lnk){ return lnk.dist; })
            .method("isequal", &Link::operator==)
            .method("isless", &Link::operator<)
            .method("<<", [](Link& l){ std::cout << l; });

    mod.add_type<jlcxx::Parametric<jlcxx::TypeVar<1>>>("BSGVector")
        .apply<std::vector<Node>,
                std::vector<sgNodeID_t>,
                std::vector<Link>,
                std::vector<nodeVisitor>,
                std::vector<SequenceGraphPath>,
                std::vector<SequenceSubGraph>,
                std::vector<std::vector<sgNodeID_t>>>([](auto wrapped){
            typedef typename decltype(wrapped)::type WrappedT;
            wrapped.method("size", &WrappedT::size);
            wrapped.method("at", [](WrappedT& vec, size_t i){ return vec.at(i); });
    });
/*
    mod.add_type<NodeIDVec>("NodeIDVec")
            .method("size", &NodeIDVec::size)
            .method("at", [](NodeIDVec& vec, size_t i){
                return vec.at(i);
            });
    mod.add_type<LinkVec>("LinkVec")
            .method("size", &LinkVec::size)
            .method("at", [](LinkVec& vec, size_t i){
               return vec.at(i);
            });
    mod.add_type<NodeVec>("NodeVec")
            .method("size", &NodeVec::size)
            .method("at", [](NodeVec& vec, size_t i){
                return vec.at(i);
            });
            */

    mod.add_type<jlcxx::Parametric<jlcxx::TypeVar<1>>>("BSGSet")
        .apply<std::set<nodeVisitor>>([](auto wrapped){
            typedef typename decltype(wrapped)::type WrappedT;
            wrapped.constructor();
            wrapped.method("size", &WrappedT::size);
            wrapped.method("clear", &WrappedT::clear);
            //wrapped.method("insert", &WrappedT::insert);
            //wrapped.method("erase", &WrappedT::erase);
        });

    // Graph types
    mod.add_type<SequenceGraph>("SequenceGraph")
            .constructor()
            .method("load_from_gfa", &SequenceGraph::load_from_gfa)
            .method("write_to_gfa", [](SequenceGraph& sg, std::string filename) -> void {
                sg.write_to_gfa(filename);
            })
            .method("add_node", &SequenceGraph::add_node)
            .method("get_node", &SequenceGraph::get_node)
            .method("add_link", &SequenceGraph::add_link)
            .method("get_link", &SequenceGraph::get_link)
            .method("get_fw_links", &SequenceGraph::get_fw_links)
            .method("get_bw_links", &SequenceGraph::get_bw_links)
            .method("get_fw_nodes", &SequenceGraph::get_fw_nodes)
            .method("get_bw_nodes", &SequenceGraph::get_bw_nodes)
            .method("get_neighbour_nodes", &SequenceGraph::get_neighbour_nodes)
            .method("is_sane", &SequenceGraph::is_sane)
            .method("get_loopy_nodes", &SequenceGraph::get_loopy_nodes)
            .method("get_flanking_nodes", &SequenceGraph::get_flanking_nodes)
            .method("depth_first_search", &SequenceGraph::depth_first_search)
            .method("breath_first_search", &SequenceGraph::breath_first_search)
            .method("find_path_between", &SequenceGraph::find_path_between)
            .method("explore_nodes", &SequenceGraph::explore_nodes)
            .method("remove_node", &SequenceGraph::remove_node)
            .method("remove_link", &SequenceGraph::remove_link)
            .method("join_path", &SequenceGraph::join_path)
            .method("join_all_unitigs", &SequenceGraph::join_all_unitigs)
            //.method("get_all_unitigs", &SequenceGraph::get_all_unitigs)
            //.method("expand_node", &SequenceGraph::expand_node)
                    //.method("get_distances_to", &SequenceGraph::get_distances_to)
            //.method("get_all_bubbly_subgraphs", &SequenceGraph::get_all_bubbly_subgraphs)
            //.method("print_bubbly_subgraph_stats", &SequenceGraph::print_bubbly_subgraph_stats)
            .method("oldnames_to_nodes", &SequenceGraph::oldnames_to_nodes)
            .method("get_node", &SequenceGraph::get_node)
                    //.method("is_loop", &SequenceGraph::is_loop)
            .method("link_exists", &SequenceGraph::link_exists)
            .method("nodeID_to_name", &SequenceGraph::nodeID_to_name)
            .method("count_active_nodes", &SequenceGraph::count_active_nodes);
            //.method("find_all_paths_between", &SequenceGraph::find_all_paths_between);

    mod.add_type<SequenceGraphPath>("SequenceGraphPath")
            .constructor<SequenceGraph&, const std::vector<sgNodeID_t>>()
            .constructor<const SequenceGraph&, const std::vector<sgNodeID_t>>()
            .constructor<const SequenceGraphPath&>()
            .method("get_fasta_header", &SequenceGraphPath::get_fasta_header)
            .method("get_sequence", &SequenceGraphPath::get_sequence)
            .method("get_sequence_size_fast", &SequenceGraphPath::get_sequence_size_fast)
            .method("get_next_links", &SequenceGraphPath::get_next_links)
            .method("reverse", &SequenceGraphPath::reverse)
            .method("is_canonical", &SequenceGraphPath::is_canonical)
            //.method("make_set_of_nodes", &SequenceGraphPath::make_set_of_nodes)
            .method("append_to_path", &SequenceGraphPath::append_to_path)
            .method("extend_if_coherent", &SequenceGraphPath::extend_if_coherent)
            .method("clear", &SequenceGraphPath::clear)
            .method("is_unitig", &SequenceGraphPath::is_unitig)
            .method("getNodes", static_cast<std::vector<sgNodeID_t>& (SequenceGraphPath::*)()>(&SequenceGraphPath::getNodes));

 /*
    mod.add_type<SequenceSubGraph>("SequenceSubGraph")
            .constructor<const SequenceGraph&, std::vector<sgNodeID_t>>()
            .constructor<const SequenceGraph&, std::vector<nodeVisitor>>()
            .method("write_to_gfa", &SequenceSubGraph::write_to_gfa)
            .method("total_size", &SequenceSubGraph::total_size);
*/
 /*
    mod.add_type<LinkageDiGraph>("LinkageDiGraph")
            .constructor<SequenceGraph&>()
            .method("add_link", &LinkageDiGraph::add_link)
            .method("add_links", &LinkageDiGraph::add_links)
            .method("remove_link", &LinkageDiGraph::remove_link)
            .method("disconnect_node", &LinkageDiGraph::disconnect_node)
            .method("get_fw_links", &LinkageDiGraph::get_fw_links)
            .method("get_bw_links", &LinkageDiGraph::get_bw_links)
            .method("fw_reached_nodes", &LinkageDiGraph::fw_reached_nodes)
            .method("get_connected_nodes", &LinkageDiGraph::get_connected_nodes)
            .method("remove_transitive_links", &LinkageDiGraph::remove_transitive_links)
            .method("report_connectivity", &LinkageDiGraph::report_connectivity)
            .method("are_connected", &LinkageDiGraph::are_connected)
                    //.method("get_all_lines", &LinkageDiGraph::get_all_lines)
            .method("find_bubbles", &LinkageDiGraph::find_bubbles)
            .method("dump_to_text", &LinkageDiGraph::dump_to_text)
            .method("load_from_text", &LinkageDiGraph::load_from_text);
*/



/*
    // Untanglers
    mod.add_type<LinkageUntangler>("LinkageUntangler")
            .constructor<WorkSpace&>()
            .method("clear_node_selection", &LinkageUntangler::clear_node_selection)
            .method("report_node_selection", &LinkageUntangler::report_node_selection)
            .method("select_nodes_by_size_and_ci", &LinkageUntangler::select_nodes_by_size_and_ci)
            .method("get_HSPNPs", &LinkageUntangler::get_HSPNPs)
            .method("select_nodes_by_HSPNPs", &LinkageUntangler::select_nodes_by_HSPNPs)
            .method("shared_read_paths", &LinkageUntangler::shared_read_paths)
            .method("make_topology_linkage", &LinkageUntangler::make_topology_linkage)
            .method("make_paired_linkage", &LinkageUntangler::make_paired_linkage)
            .method("make_longRead_linkage", &LinkageUntangler::make_longRead_multilinkage)
            .method("make_paired_linkage_pe", &LinkageUntangler::make_paired_linkage_pe)
            .method("make_paired_linkage_by_kmer", &LinkageUntangler::make_paired_linkage_by_kmer)
            .method("make_tag_linkage", &LinkageUntangler::make_tag_linkage)
            .method("make_and_simplify_linkage", &LinkageUntangler::make_and_simplify_linkage)
            .method("filter_linkage_to_hspnp_duos", &LinkageUntangler::filter_linkage_to_hspnp_duos)
            .method("expand_trivial_repeats", &LinkageUntangler::expand_trivial_repeats)
            .method("expand_linear_regions", &LinkageUntangler::expand_linear_regions)
            .method("expand_linear_regions_skating", &LinkageUntangler::expand_linear_regions_skating)
            .method("linear_regions_tag_local_assembly", &LinkageUntangler::linear_regions_tag_local_assembly)
            .method("fill_linkage_line", &LinkageUntangler::fill_linkage_line);
            */

    /*
    // Read mappers
    mod.add_type<LinkedReadMapper>("LinkedReadMapper")
            .constructor<SequenceGraph&, LinkedReadsDatastore&, const UniqueKmerIndex&, const Unique63merIndex&>()
            .method("map_reads", &LinkedReadMapper::map_reads)
            .method("remap_all_reads", &LinkedReadMapper::remap_all_reads)
            .method("map_reads63", &LinkedReadMapper::map_reads63)
            .method("remap_all_reads63", &LinkedReadMapper::remap_all_reads63)
            .method("remove_obsolete_mappings", &LinkedReadMapper::remove_obsolete_mappings)
            .method("get_node_tags", &LinkedReadMapper::get_node_tags)
            .method("get_tag_nodes", &LinkedReadMapper::get_tag_nodes)
            .method("get_tag_neighbour_nodes", &LinkedReadMapper::get_tag_neighbour_nodes)
            .method("print_stats", &LinkedReadMapper::print_stats); */
/*
    mod.add_type<PairedReadMapper>("PairedReadMapper")
            .constructor<SequenceGraph&, PairedReadsDatastore&, const UniqueKmerIndex&, const Unique63merIndex&>()
            .method("write", &PairedReadMapper::write)
            .method("read", &PairedReadMapper::read)
            .method("map_reads", &PairedReadMapper::map_reads)
            .method("remap_all_reads", &PairedReadMapper::remap_all_reads)
            .method("map_reads63", &PairedReadMapper::map_reads63)
            .method("remap_all_reads63", &PairedReadMapper::remap_all_reads63)
            .method("remove_obsolete_mappings", &PairedReadMapper::remove_obsolete_mappings)
            .method("size_distribution", &PairedReadMapper::size_distribution)
            .method("populate_orientation", &PairedReadMapper::populate_orientation)
            .method("print_stats", &PairedReadMapper::print_stats)
            .method("get_node_readpairs_ids", &PairedReadMapper::get_node_readpairs_ids);

            */
/*
    mod.add_type<LongReadMapper>("LongReadMapper")
            .constructor<SequenceGraph&, LongReadsDatastore&, uint8_t>()
            .method("getLongReadsDatastore", &LongReadMapper::getLongReadsDatastore)
            .method("get_node_read_ids", &LongReadMapper::get_node_read_ids)
            .method("set_params", &LongReadMapper::set_params)
            .method("get_all_kmer_matches", &LongReadMapper::get_all_kmer_matches)
            .method("window_candidates", &LongReadMapper::window_candidates)
            .method("alignment_blocks", &LongReadMapper::alignment_blocks)
            .method("filter_blocks", &LongReadMapper::filter_blocks)
            .method("refine_multinode_reads", &LongReadMapper::refine_multinode_reads)
            .method("map_reads", static_cast<void (LongReadMapper::*)(std::string)>(&LongReadMapper::map_reads))
            .method("map_reads", static_cast<void (LongReadMapper::*)(std::unordered_set<uint32_t>, std::string)>(&LongReadMapper::map_reads))
            .method("read", static_cast<void (LongReadMapper::*)(std::string)>(&LongReadMapper::read))
            .method("read", static_cast<void (LongReadMapper::*)(std::ifstream&)>(&LongReadMapper::read))
            .method("write", static_cast<void (LongReadMapper::*)(std::string)>(&LongReadMapper::write))
            .method("write", static_cast<void (LongReadMapper::*)(std::ofstream&)>(&LongReadMapper::write))
            .method("update_graph_index", &LongReadMapper::update_graph_index);
*/
/*
    mod.add_type<KmerCompressionIndex>("KmerCompressionIndex")
            .constructor<SequenceGraph&, uint64_t>()
            .constructor<SequenceGraph&>()
            .method("index_graph", &KmerCompressionIndex::index_graph)
            .method("reindex_graph", &KmerCompressionIndex::reindex_graph)
            .method("start_new_count", &KmerCompressionIndex::start_new_count)
            .method("add_counts_from_file", &KmerCompressionIndex::add_counts_from_file)
            .method("add_counts_from_datastore", &KmerCompressionIndex::add_counts_from_datastore)
            .method("write", &KmerCompressionIndex::write)
            .method("read", &KmerCompressionIndex::read)
            .method("save_to_disk", &KmerCompressionIndex::save_to_disk)
            .method("load_from_disk", &KmerCompressionIndex::load_from_disk)
            .method("compute_compression_stats", &KmerCompressionIndex::compute_compression_stats)
            //.method("compute_all_nodes_kci", &KmerCompressionIndex::compute_all_nodes_kci)
            .method("compute_kci_profiles", &KmerCompressionIndex::compute_kci_profiles)
            .method("compute_node_coverage_profile", &KmerCompressionIndex::compute_node_coverage_profile)
            .method("dump_histogram", &KmerCompressionIndex::dump_histogram);
            //.method("compute_compression_for_node", &KmerCompressionIndex::compute_compression_for_node);
*/

    // Workspaces
    mod.add_type<WorkSpace>("WorkSpace")
            .constructor()
            .method("print_log", &WorkSpace::print_log)
            .method("add_log_entry", &WorkSpace::add_log_entry)
            .method("load_from_disk", &WorkSpace::load_from_disk)
            .method("dump_to_disk", &WorkSpace::dump_to_disk)
            //.method("getKCI", &WorkSpace::getKCI)
            .method("getGraph", &WorkSpace::getGraph);
            //.method("getLinkedReadMappers", &WorkSpace::getLinkedReadMappers)
            //.method("getLinkedReadDatastores", &WorkSpace::getLinkedReadDatastores)
            //.method("getPairedReadMappers", &WorkSpace::getPairedReadMappers)
            //.method("getPairedReadDatastores", &WorkSpace::getPairedReadDatastores)
            //.method("getLongReadMappers", &WorkSpace::getLongReadMappers);
            //.method("getLongReadDatastores", &WorkSpace::getLongReadDatastores)
            //.method("getPathsDatastore", &WorkSpace::getPathsDatastore);

}