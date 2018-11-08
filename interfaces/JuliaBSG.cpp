//
// Created by Ben Ward (EI) on 02/10/2018.
//

#include "jlcxx/jlcxx.hpp"
#include <sglib/graph/SequenceGraph.hpp>
#include <sglib/workspace/WorkSpace.hpp>
#include <sglib/processors/LinkageUntangler.hpp>
#include <sglib/mappers/PairedReadMapper.hpp>
#include <sglib/graph/LinkageDiGraph.hpp>

JLCXX_MODULE define_julia_module(jlcxx::Module& mod)
{
    mod.add_type<SequenceGraph>("SequenceGraph")
            .constructor()
            .method("load_from_gfa", &SequenceGraph::load_from_gfa)
            .method("write_to_gfa", &SequenceGraph::write_to_gfa)
            .method("add_node", &SequenceGraph::add_node)
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
            .method("oldnames_to_nodes", &SequenceGraph::oldnames_to_nodes)
            .method("get_node", &SequenceGraph::get_node)
            .method("link_exists", &SequenceGraph::link_exists)
            .method("nodeID_to_name", &SequenceGraph::nodeID_to_name);

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
            .method("make_longRead_linkage", &LinkageUntangler::make_longRead_linkage)
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
            .method("print_stats", &LinkedReadMapper::print_stats);

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

    mod.add_type<WorkSpace>("WorkSpace")
            .constructor()
            .method("load_from_disk", &WorkSpace::load_from_disk)
            .method("dump_to_disk", &WorkSpace::dump_to_disk)
            .method("getKCI", &WorkSpace::getKCI)
            .method("getGraph", &WorkSpace::getGraph)
            .method("getLinkedReadMappers", &WorkSpace::getLinkedReadMappers)
            .method("getLinkedReadDatastores", &WorkSpace::getLinkedReadDatastores)
            .method("getPairedReadMappers", &WorkSpace::getPairedReadMappers)
            .method("getPairedReadDatastores", &WorkSpace::getPairedReadDatastores)
            .method("getLongReadMappers", &WorkSpace::getLongReadMappers)
            .method("getLongReadDatastores", &WorkSpace::getLongReadDatastores)
            .method("getPathsDatastore", &WorkSpace::getPathsDatastore);

}