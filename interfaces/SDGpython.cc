//
// Created by Bernardo Clavijo (EI) on 13/02/2020.
//

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>
#include <sdglib/graph/SequenceDistanceGraph.hpp>
#include <sdglib/workspace/WorkSpace.hpp>
#include <sdglib/views/NodeView.hpp>
#include <sdglib/views/TangleView.hpp>
#include <sdglib/mappers/LongReadsRecruiter.hpp>
#include <sdglib/processors/GraphEditor.hpp>
#include <sdglib/processors/LinkageUntangler.hpp>
#include <sdglib/processors/GraphContigger.hpp>
#include <sdglib/processors/GraphMaker.hpp>
#include <sdglib/processors/LinkageMaker.hpp>
#include <sdglib/processors/GraphPatcher.hpp>
#include <sdglib/processors/PathFinder.hpp>
#include <sdglib/processors/Strider.hpp>
#include <sdglib/processors/CountFilter.hpp>
#include <sdglib/batch_counter/BatchKmersCounter.hpp>
#include <sdglib/processors/LineFiller.h>
#include <sdglib/processors/HappySorter.hpp>
#include <sdglib/graph/ReadThreadsGraph.hpp>
#include <sdglib/processors/TotalSorter.hpp>
#include <sdglib/processors/RTGCluster.hpp>
#include <sdglib/processors/RTGClassifier.hpp>

namespace py = pybind11;
using namespace py::literals;


PYBIND11_MAKE_OPAQUE(std::vector<std::vector<PerfectMatch>>);
PYBIND11_MAKE_OPAQUE(std::vector<std::vector<uint64_t>>);
PYBIND11_MAKE_OPAQUE(std::vector<std::vector<int64_t>>);
PYBIND11_MAKE_OPAQUE(std::vector<std::vector<NodePosition>>);
PYBIND11_MAKE_OPAQUE(std::vector<std::vector<std::pair<uint32_t,uint32_t>>>);
PYBIND11_MAKE_OPAQUE(std::vector<ReadPath>);
PYBIND11_MAKE_OPAQUE(std::vector<bool>);
PYBIND11_MAKE_OPAQUE(std::vector<std::vector<Link>>);
PYBIND11_MAKE_OPAQUE(std::vector<NodeView>);
PYBIND11_MAKE_OPAQUE(std::vector<TangleView>);

PYBIND11_MODULE(SDGpython, m) {

    py::bind_vector<std::vector<uint64_t>>(m, "VectorU64");
    py::bind_vector<std::vector<std::vector<uint64_t>>>(m, "VectorVectorU64");
    py::bind_vector<std::vector<int64_t>>(m, "Vector64");
    py::bind_vector<std::vector<std::vector<int64_t>>>(m, "VectorVector64");
    py::bind_vector<std::vector<std::vector<std::pair<uint32_t,uint32_t>>>>(m, "VectorVectorPairU32");
    py::bind_vector<std::vector<ReadPath>>(m, "VectorReadPath");
    py::bind_vector<std::vector<PerfectMatch>>(m, "VectorPerfectMatch");
    py::bind_vector<std::vector<std::vector<PerfectMatch>>>(m, "VectorVectorPerfectMatch");
    py::bind_vector<std::vector<std::vector<Link>>>(m, "VectorVectorLink");
    py::bind_vector<std::vector<std::vector<NodePosition>>>(m, "VectorVectorNodePosition");
    py::bind_vector<std::vector<bool>>(m, "VectorBool");
    py::bind_vector<std::vector<NodeView>>(m, "VectorNodeView");
    py::bind_vector<std::vector<TangleView>>(m, "VectorTangleView");
    py::bind_map<std::map<uint64_t,std::vector<std::pair<int64_t,sgNodeID_t>>>>(m,"MapThreadPosNode");


    py::enum_<NodeStatus>(m,"NodeStatus")
            .value("Active",NodeStatus::Active)
            .value("Deleted",NodeStatus::Deleted)
            .export_values();

    py::enum_<SupportType >(m,"SupportType")
            .value("stUndefined",SupportType::Undefined)
            .value("stOperation",SupportType::Operation)
            .value("stSequenceDistanceGraph",SupportType::SequenceDistanceGraph)
            .value("stDistanceGraph",SupportType::DistanceGraph)
            .value("stPairedRead",SupportType::PairedRead)
            .value("stLinkedRead",SupportType::LinkedRead)
            .value("stLinkedTag",SupportType::LinkedTag)
            .value("stLongRead",SupportType::LongRead)
            .value("stReadPath",SupportType::ReadPath)
            .value("stKmerCoverage",SupportType::KmerCoverage)
            .value("stStriderLink",SupportType::StriderLink)
            .export_values();

    py::enum_<KmerCountMode>(m,"KmerCountMode")
            .value("Canonical",KmerCountMode::Canonical)
            .value("NonCanonical",KmerCountMode::NonCanonical)
            .export_values();

    py::class_<Support>(m,"Support","Support for an operation or piece of data")
            .def(py::init<SupportType,uint16_t,int64_t>(),"","support"_a=SupportType::Undefined,"index"_a=0,"id"_a=0,py::return_value_policy::take_ownership)
            .def_readonly("type",&Support::type)
            .def_readonly("index",&Support::index)
            .def_readonly("id",&Support::id)
            ;
    py::class_<Node>(m, "Node", "A node in a Sequence Distance Graph")
            .def(py::init<const std::string &>(),"","sequence"_a,py::return_value_policy::take_ownership)
            .def_readonly("sequence", &Node::sequence)
            .def("is_canonical", &Node::is_canonical)
            .def_readonly("status", &Node::status)
            .def("__repr__",
                 [](const Node &n) {
                     return "<Node with " + std::to_string(n.sequence.size()) + " bp>";
                 })
        //.def_readonly("support", &Node::support)
            ;

    py::class_<Link>(m,"Link","A raw DG Link")
            .def(py::init<sgNodeID_t, sgNodeID_t, int32_t, Support>(),"","src"_a,"dest"_a,"dist"_a,"support"_a=Support(),py::return_value_policy::take_ownership)
            .def_readwrite("source",&Link::source)
            .def_readwrite("dest",&Link::dest)
            .def_readwrite("dist",&Link::dist)
            .def_readwrite("support",&Link::support)
            .def("__repr__",
                 [](const Link &l) {
                     return "<Link ( " + std::to_string(l.source) + " , " + std::to_string(l.dest) + " ), d = " + std::to_string(l.dist) + " bp>";
                 })
            ;

    py::class_<LinkView>(m,"LinkView","A view for a Link in a Distance Graph")
            .def("node",&LinkView::node)
            .def("distance",&LinkView::distance)
            .def("support",&LinkView::support)
            .def("__repr__",
                 [](const LinkView &lv) {
                     return "<LinkView to " + std::to_string(lv.node().node_id()) + " at " + std::to_string(lv.distance()) + " bp>";
                 })
            ;

    py::class_<NodeView>(m,"NodeView","A view for a Node in a Distance Graph")
            .def("node_id",&NodeView::node_id)
            .def("rc",&NodeView::rc)
            .def("is_tip",&NodeView::is_tip)
            .def("is_bubble_side",&NodeView::is_bubble_side)
            .def("is_canonical_repeat",&NodeView::is_canonical_repeat)
            .def("sequence",&NodeView::sequence)
            .def("size",&NodeView::size)
            .def("prev",&NodeView::prev,"A list with LinkViews to previous nodes",py::return_value_policy::move)
            .def("next",&NodeView::next,"A list with LinkViews to next nodes",py::return_value_policy::move)
            .def("parallels",&NodeView::parallels,"A list with NodeViews of parallel nodes",py::return_value_policy::move)
            .def("kmer_coverage",py::overload_cast<std::string, std::string>(&NodeView::kmer_coverage, py::const_))
            .def("kci",&NodeView::kci)
            .def("get_mappings_tags", &NodeView::get_mappings_tags,"dsname"_a,"min_reads"_a=3)
            .def("get_paths_tags", &NodeView::get_paths_tags,"dsname"_a,"min_reads"_a=3)
            .def("get_kmers",&NodeView::get_kmers, "K"_a=31)
            .def("__eq__", &NodeView::operator==, py::is_operator())
            .def("__repr__",
                 [](const NodeView &nv) {
                     return "<NodeView: node " + std::to_string(nv.node_id()) + " in graph " + nv.graph().name + ">";
                 })
            ;

    py::class_<TangleView>(m,"TangleView","A view for an unresolved section in a Distance Graph")
            .def_readonly("internals",&TangleView::internals)
            .def_readonly("frontiers",&TangleView::frontiers)
            .def("classify_tangle",&TangleView::classify_tangle)
            .def("__repr__",
                 [](const TangleView &tv) {
                     return "<TangleView: " + std::to_string(tv.internals.size()) + " internals and " + std::to_string(tv.frontiers.size()) + " frontiers in graph " + tv.dg->name + ">";
                 })
            ;

    py::class_<DistanceGraph>(m, "DistanceGraph", "A Distance Graph")
            .def_property_readonly("sdg", [] (const DistanceGraph &dg) {return &dg.sdg;},py::return_value_policy::reference)
            .def(py::init<SequenceDistanceGraph &, const std::string &>(),"","sdg"_a,"name"_a="unnamed",py::return_value_policy::take_ownership)
            .def("add_link",&DistanceGraph::add_link,"source"_a,"dest"_a,"distance"_a,"support"_a=Support() )
            .def("disconnect_node",&DistanceGraph::disconnect_node)
            .def("remove_link",py::overload_cast<sgNodeID_t , sgNodeID_t >(&DistanceGraph::remove_link))
            .def("remove_link",py::overload_cast<sgNodeID_t , sgNodeID_t, int32_t, Support >(&DistanceGraph::remove_link))
            .def("remove_transitive_links",&DistanceGraph::remove_transitive_links,"radius"_a=10)
            .def("fw_neighbours_by_distance",&DistanceGraph::fw_neighbours_by_distance,"node_id"_a,"min_links"_a=3)
            .def("fw_reached_nodes",&DistanceGraph::fw_reached_nodes,"node_id"_a,"radius"_a=10)
            .def("get_nodeview",&DistanceGraph::get_nodeview)
            .def("get_all_nodeviews",&DistanceGraph::get_all_nodeviews,"both_directions"_a=false,"include_disconnected"_a=true,"Returns a vector with NodeViews for all active nodes",py::return_value_policy::move)
            .def("get_all_tangleviews",&DistanceGraph::get_all_tangleviews,"f_size"_a=500,"f_min_kci"_a=-1,"f_max_kci"_a=-1,"include_disconnected"_a=true,"Returns a vector with TangleViews for all tangles",py::return_value_policy::move)
            .def("get_all_anchors_tangleviews",&DistanceGraph::get_all_anchors_tangleviews, "given_frontiers"_a, "include_disconnected"_a=true,"Returns a vector with TangleViews for all tangles",py::return_value_policy::move)
            .def_readwrite("name",&DistanceGraph::name)
            .def("write_to_gfa1",&DistanceGraph::write_to_gfa1,"filename"_a,"selected_nodes"_a=std::vector<sgNodeID_t>(),"depths"_a=std::vector<sgNodeID_t>())
            .def("write_to_gfa2",&DistanceGraph::write_to_gfa2)
            .def("stats_by_kci",&DistanceGraph::stats_by_kci)
            .def("simple_structure_stats",&DistanceGraph::simple_structure_stats)
            .def("are_connected",&DistanceGraph::are_connected,"node1"_a, "node2"_a)
            .def("find_all_paths_between",&DistanceGraph::find_all_paths_between,"from"_a,"to"_a,"max_size"_a,"max_nodes"_a,"abort_on_loops"_a)
            .def("get_all_lines",&DistanceGraph::get_all_lines,"min_nodes"_a, "min_total_size"_a=0)
            .def("get_connected_component",&DistanceGraph::get_connected_component,"node"_a, "signed_nodes"_a=true)
            .def("get_all_connected_components",&DistanceGraph::get_all_connected_components,"min_nodes"_a=1)
            ;

    py::class_<SequenceDistanceGraph,DistanceGraph>(m, "SequenceDistanceGraph", "A Sequence Distance Graph")
            .def("get_node_size",&SequenceDistanceGraph::get_node_size)
            .def("add_node",py::overload_cast<std::string,bool>(&SequenceDistanceGraph::add_node),"sequence"_a,"make_canonical"_a=true)
            .def("remove_node",&SequenceDistanceGraph::remove_node)
            .def("join_all_unitigs",&SequenceDistanceGraph::join_all_unitigs)
            .def("load_from_gfa",&SequenceDistanceGraph::load_from_gfa)
            .def("load_from_bcalm",&SequenceDistanceGraph::load_from_bcalm,"filename"_a,"k"_a)
            .def("load_from_fasta",&SequenceDistanceGraph::load_from_fasta)
            ;

    py::class_<ReadThreadsGraph,DistanceGraph>(m, "ReadThreadsGraph", "A Read Threads Graph")
            .def(py::init<SequenceDistanceGraph &, const std::string &>(),"","sdg"_a,"name"_a="unnamed",py::return_value_policy::take_ownership)
            .def("dump",&ReadThreadsGraph::dump, "filename"_a)
            .def("load",&ReadThreadsGraph::load, "filename"_a)
            .def("add_thread",&ReadThreadsGraph::add_thread, "thread_id"_a, "thread"_a,"remove_duplicated"_a=true,"min_nodes"_a=2)
            .def("remove_thread",&ReadThreadsGraph::remove_thread, "thread_id"_a)
            .def("get_thread",&ReadThreadsGraph::get_thread, "thread_id"_a)
            .def("thread_start_nodeview", &ReadThreadsGraph::thread_start_nodeview, "thread_id"_a)
            .def("thread_end_nodeview", &ReadThreadsGraph::thread_end_nodeview, "thread_id"_a)
            .def("next_in_thread",&ReadThreadsGraph::next_in_thread, "nid"_a,"thread_id"_a,"link_index"_a=-1)
            .def("prev_in_thread",&ReadThreadsGraph::prev_in_thread, "nid"_a,"thread_id"_a,"link_index"_a=-1)
            .def("all_nids_fw_in_thread",&ReadThreadsGraph::all_nids_fw_in_thread, "nid"_a,"thread_id"_a)
            .def("local_graph",&ReadThreadsGraph::local_graph, "node_id"_a,"distance"_a,"min_links"_a)
            .def("pop_node",&ReadThreadsGraph::pop_node,"node_id"_a,"thread_id"_a)
            .def("pop_nodes",&ReadThreadsGraph::pop_nodes,"node_ids"_a,"thread_id"_a)
            .def("pop_node_from_all",&ReadThreadsGraph::pop_node_from_all,"node_id"_a)
            .def("flip_thread",&ReadThreadsGraph::flip_thread,"thread_id"_a)
            .def("node_threads",&ReadThreadsGraph::node_threads,"node_id"_a,"oriented"_a=false)
            .def("node_threadpositions",&ReadThreadsGraph::node_threadpositions,"node_id"_a)
            .def("node_thread_neighbours",&ReadThreadsGraph::node_thread_neighbours,"node_id"_a,"oriented"_a=false)
            .def("clean_node",&ReadThreadsGraph::clean_node,"node_id"_a,"min_supported"_a=4,"min_support"_a=1)
            .def("clean_lowmapping_nodes_popping_list",&ReadThreadsGraph::clean_lowmapping_nodes_popping_list,"min_threads"_a=4)
            .def("clean_repeat_nodes_popping_list",&ReadThreadsGraph::clean_repeat_nodes_popping_list,"max_threads"_a=200)
            .def("clean_all_nodes_popping_list",&ReadThreadsGraph::clean_all_nodes_popping_list,"min_supported"_a=4,"min_support"_a=1)
            .def("clean_all_nodes_by_thread_clustering_popping_list",&ReadThreadsGraph::clean_all_nodes_by_thread_clustering_popping_list,"min_shared"_a=4,"max_second_perc"_a=.1)
            .def("thread_nodesets",&ReadThreadsGraph::thread_nodesets)
            .def("apply_popping_list",&ReadThreadsGraph::apply_popping_list,"popping_list"_a)
            .def("thread_fw_in_node",&ReadThreadsGraph::thread_fw_in_node,"thread_id"_a,"node_id"_a)
            .def("nodes_before_in_thread",&ReadThreadsGraph::nodes_before_in_thread,"thread_id"_a,"node_id"_a)
            .def("nodes_after_in_thread",&ReadThreadsGraph::nodes_after_in_thread,"thread_id"_a,"node_id"_a)
            .def("make_thread_nodepositions",&ReadThreadsGraph::make_thread_nodepositions,"nodes"_a)
            ;


    py::class_<HappySorter>(m, "HappySorter", "HappySorter")
            .def(py::init<const ReadThreadsGraph &, float, int, int>(),"","rtg"_a, "min_node_happiness"_a=.7,
             "min_node_threads"_a=2, "order_end_size"_a=20,py::return_value_policy::take_ownership)
            .def("node_happiness", &HappySorter::node_happiness,"tid"_a,"prev"_a=true,"next"_a=false,"min_threads"_a=3)
            .def("start_from_nodelist",&HappySorter::start_from_nodelist,"nodelist"_a,"min_links"_a=3)
            .def("find_fw_candidates",&HappySorter::find_fw_candidates,"min_happiness"_a=-1,"min_threads"_a=-1,"end_size"_a=-1)
            .def("find_bw_candidates",&HappySorter::find_bw_candidates,"min_happiness"_a=-1,"min_threads"_a=-1,"end_size"_a=-1)
            .def("find_internal_candidates",&HappySorter::find_internal_candidates,"min_happiness"_a=-1,"min_threads"_a=-1,"first"_a=1,"last"_a=INT32_MAX)
            .def("close_internal_threads",&HappySorter::close_internal_threads,"order_end"_a=30,"thread_end"_a=0)
            .def("make_thread_nodepositions",&HappySorter::make_thread_nodepositions,"nodes"_a,"tids"_a=std::set<int64_t>())
            .def("place_nodes",&HappySorter::place_nodes,"nodes"_a,"verbose"_a=false)
            .def("thread_happiness_q",&HappySorter::thread_happiness_q,"tid"_a,"min_nodes"_a,"max_span"_a)
            .def("recruit_all_happy_threads_q", &HappySorter::recruit_all_happy_threads_q,"min_nodes"_a=7,"max_span"_a=10)
            .def("q_grow_loop",&HappySorter::q_grow_loop, "min_threads"_a=-1, "min_happiness"_a=-1, "p"_a=4, "q"_a=5, "steps"_a=INT64_MAX, "verbose"_a=false)
            .def("grow",&HappySorter::grow, "_thread_hits"_a=4, "_end_size"_a=30, "_node_hits"_a=3, "_min_happiness"_a=.7, "verbose"_a=false)
            .def("run_from_nodelist",&HappySorter::run_from_nodelist, "nodes"_a,"min_threads"_a=-1, "min_happiness"_a=-1, "p"_a=4, "q"_a=5, "steps"_a=INT64_MAX, "verbose"_a=false)
            .def("hs_tnp_to_distances",&HappySorter::hs_tnp_to_distances, "thread_nodepositions"_a, "nodeset"_a)
            .def("update_positions",&HappySorter::update_positions,"first"_a=0,"last"_a=-1)
            .def("is_mixed",&HappySorter::is_mixed,"win"_a=50,"fail"_a=.2)
            .def_readwrite("threads",&HappySorter::threads)
            .def_readwrite("order",&HappySorter::order)
            .def_readwrite("bw_open_threads",&HappySorter::bw_open_threads)
            .def_readwrite("fw_open_threads",&HappySorter::fw_open_threads)
            ;

    pybind11::class_<TotalSorter>(m, "TotalSorter", "TotalSorter")
            .def(py::init<ReadThreadsGraph &, int, int>(),"","rtg"_a,"min_thread_length"_a=6, "min_node_threads"_a=3,py::return_value_policy::take_ownership)
            .def(py::init<ReadThreadsGraph &, std::string>(),"","rtg"_a,"filename"_a,py::return_value_policy::take_ownership)
            .def("run_sorters_from_lines",&TotalSorter::run_sorters_from_lines,"lines"_a=std::vector<std::vector<sgNodeID_t>>(),"min_line_size"_a=25,"min_order_size"_a=200,"line_occupancy"_a=.7, "p"_a=5, "q"_a=10, "node_happiness"_a=.8, "node_threads"_a=5, "end_size"_a=30,"min_links"_a=5)
            .def("prune_rtg",&TotalSorter::prune_rtg)
            .def("update_usage",&TotalSorter::update_usage)
            .def("merge",&TotalSorter::merge)
            .def("remove_mixed",&TotalSorter::remove_mixed,"win"_a=50,"fail"_a=.2)
            .def("compute_node_neighbours",&TotalSorter::compute_node_neighbours,"k"_a,"max_f"_a)
            .def("dump",&TotalSorter::dump,"filename"_a)
            .def("load",&TotalSorter::load,"filename"_a)
            .def_readwrite("nodes",&TotalSorter::nodes)
            .def_readwrite("threads",&TotalSorter::threads)
            .def_readonly("sorters",&TotalSorter::sorters)
            .def_readonly("merged_sorters",&TotalSorter::sorters)
            .def_readonly("sorter_classes",&TotalSorter::sorter_classes)
            .def_readonly("node_sorters",&TotalSorter::node_sorters)
            .def_readonly("thread_sorters",&TotalSorter::thread_sorters)
            .def_readonly("node_neighbours",&TotalSorter::node_neighbours)
                 ;

    py::class_<SequenceDistanceGraphPath>(m, "SequenceDistanceGraphPath", "SequenceDistanceGraphPath")
            .def(py::init<const SequenceDistanceGraph &,const std::vector<sgNodeID_t> >(),"","sdg"_a,"nodes"_a,py::return_value_policy::take_ownership)
            .def("sequence", &SequenceDistanceGraphPath::sequence)
            .def_readwrite("nodes", &SequenceDistanceGraphPath::nodes)
            .def("__repr__",
                 [](const SequenceDistanceGraphPath &sdgp) {
                     return "<SDGPath: length: " + std::to_string(sdgp.nodes.size()) + " nodes>";
                 })
            ;

    py::class_<PairedReadsDatastore>(m, "PairedReadsDatastore", "A Paired Reads Datastore")
            .def("size",&PairedReadsDatastore::size)
            .def("get_read_sequence",&PairedReadsDatastore::get_read_sequence,"read_id"_a)
            .def("get_read_pair",&PairedReadsDatastore::get_read_pair,"read_id"_a)
            .def_readwrite("mapper",&PairedReadsDatastore::mapper)
            ;

    py::class_<ReadPath>(m, "ReadPath", "A Paired Reads ReadPath")
            .def_readonly("offset",&ReadPath::offset)
            .def_readonly("path",&ReadPath::path)
            ;

    py::class_<PairedReadsMapper>(m, "PairedReadsMapper", "A Paired Reads Mapper")
            .def("path_reads",&PairedReadsMapper::path_reads,"k"_a=63,"max_freq"_a=200,"fill_offsets"_a=true)
            .def_readonly("paths_in_node",&PairedReadsMapper::paths_in_node)
            .def_readonly("read_paths",&PairedReadsMapper::read_paths)
            .def_readonly("read_path_offsets",&PairedReadsMapper::read_path_offsets)
            .def("path_fw",&PairedReadsMapper::path_fw,"read_id"_a,"node"_a,"use_pair"_a=true,"collapse_pair"_a=true)
            .def("all_paths_fw",&PairedReadsMapper::all_paths_fw,"node"_a,"use_pair"_a=true,"collapse_pair"_a=true)
            .def("dump_readpaths",&PairedReadsMapper::dump_readpaths)
            .def("load_readpaths",&PairedReadsMapper::load_readpaths)
            .def("get_node_inmediate_neighbours",&PairedReadsMapper::get_node_inmediate_neighbours, "node"_a)
            .def("get_paths_in_node",&PairedReadsMapper::get_paths_in_node, "nid"_a)
            ;

    py::class_<LinkedReadsMapper>(m, "LinkedReadsMapper", "A Paired Reads Mapper")
            .def("path_reads",&LinkedReadsMapper::path_reads,"k"_a=63,"max_freq"_a=200)
            .def_readonly("paths_in_node",&LinkedReadsMapper::paths_in_node)
            .def_readonly("read_paths",&LinkedReadsMapper::read_paths)
            .def("dump_readpaths",&LinkedReadsMapper::dump_readpaths)
            .def("load_readpaths",&LinkedReadsMapper::load_readpaths)
            ;

    py::class_<LinkedReadsDatastore>(m, "LinkedReadsDatastore", "A Linked Reads Datastore")
            .def("size",&LinkedReadsDatastore::size)
            .def("get_read_sequence",&LinkedReadsDatastore::get_read_sequence)
            .def("get_read_tag",&LinkedReadsDatastore::get_read_tag)
            .def("get_read_pair",&LinkedReadsDatastore::get_read_pair,"read_id"_a)
            .def_readwrite("mapper",&LinkedReadsDatastore::mapper)
            ;

    py::class_<LongReadsDatastore>(m, "LongReadsDatastore", "A Long Reads Datastore")
            .def("size",&LongReadsDatastore::size)
            .def("get_read_sequence",&LongReadsDatastore::get_read_sequence)
            .def("get_read_size",[](const LongReadsDatastore &ds,__uint64_t rid) {
                return ds.read_to_fileRecord[rid].record_size;
            })
            .def_readwrite("mapper",&LongReadsDatastore::mapper)
            ;

    py::class_<LongReadsMapper>(m, "LongReadsMapper", "A Long Reads Mapper")
            .def("print_status",&LongReadsMapper::print_status)
            .def("set_params",&LongReadsMapper::set_params, "_k"_a=15, "_max_index_freq"_a=200, "_min_size"_a=1000, "_min_chain"_a=50, "_max_jump"_a=500, "_max_delta_change"_a= 60)
            .def("map_reads",&LongReadsMapper::map_reads, "readIDs"_a=std::unordered_set<uint32_t>())
            .def("get_reads_in_node",&LongReadsMapper::get_reads_in_node, "nid"_a)
            .def_readonly("reads_in_node",&LongReadsMapper::reads_in_node)
            .def_readonly("mappings",&LongReadsMapper::mappings)
            .def_readonly("read_paths",&LongReadsMapper::read_paths)
            .def_readwrite("k",&LongReadsMapper::k)
            ;

    py::class_<KmerCounter>(m, "KmerCounters", "A kmer counter container")
            .def("list_names",&KmerCounter::list_names)
            .def("project_count",py::overload_cast<const std::string &, const std::string & >(&KmerCounter::project_count))
            .def("set_kci_peak",&KmerCounter::set_kci_peak)
            .def("kci",&KmerCounter::kci,"node"_a)
            .def("add_count",py::overload_cast<const std::string &, const std::vector<std::string> &,bool>(&KmerCounter::add_count),"name"_a,"input_files"_a,"fastq"_a=true)
            .def("add_count",py::overload_cast<const std::string &, const PairedReadsDatastore &>(&KmerCounter::add_count),"name"_a,"datastore"_a)
            .def("add_count",py::overload_cast<const std::string &, const LinkedReadsDatastore &>(&KmerCounter::add_count),"name"_a,"datastore"_a)
            .def("add_count",py::overload_cast<const std::string &, const LongReadsDatastore &>(&KmerCounter::add_count),"name"_a,"datastore"_a)
            .def("count_spectra",&KmerCounter::count_spectra,"name"_a,"max_freq"_a=1000,"unique_in_graph"_a=false,"present_in_graph"_a=true)
            .def("update_graph_counts",&KmerCounter::update_graph_counts)
            .def("compute_all_kcis",&KmerCounter::compute_all_kcis)
            .def("dump_cache",&KmerCounter::dump_cache,"filename"_a)
            .def("load_cache",&KmerCounter::load_cache,"filename"_a)
            ;

    py::class_<WorkSpace>(m, "WorkSpace", "A full SDG WorkSpace")
            .def(py::init<const std::string &>(),"filename"_a="",py::return_value_policy::take_ownership)
            .def_readonly("sdg",&WorkSpace::sdg)
            .def("add_long_reads_datastore",&WorkSpace::add_long_reads_datastore,"datastore"_a,"name"_a="",py::return_value_policy::reference)
            .def("add_paired_reads_datastore",&WorkSpace::add_paired_reads_datastore,"datastore"_a,"name"_a="",py::return_value_policy::reference)
            .def("add_linked_reads_datastore",&WorkSpace::add_linked_reads_datastore,"datastore"_a,"name"_a="",py::return_value_policy::reference)
            .def("add_kmer_counter",py::overload_cast<const std::string &, uint8_t , KmerCountMode >(&WorkSpace::add_kmer_counter),"name"_a,"k"_a=31,"count_mode"_a=KmerCountMode::Canonical,py::return_value_policy::reference)
            .def("add_kmer_counter",py::overload_cast<const std::string &, const std::string &>(&WorkSpace::add_kmer_counter),"filename"_a,"name"_a="")
            .def("get_long_reads_datastore",&WorkSpace::get_long_reads_datastore,"name"_a,py::return_value_policy::reference)
            .def("get_paired_reads_datastore",&WorkSpace::get_paired_reads_datastore,"name"_a,py::return_value_policy::reference)
            .def("get_linked_reads_datastore",&WorkSpace::get_linked_reads_datastore,"name"_a,py::return_value_policy::reference)
            .def("get_kmer_counter",&WorkSpace::get_kmer_counter,"name"_a,py::return_value_policy::reference)
            .def("list_long_reads_datastores",&WorkSpace::list_long_reads_datastores)
            .def("list_linked_reads_datastores",&WorkSpace::list_linked_reads_datastores)
            .def("list_paired_reads_datastores",&WorkSpace::list_paired_reads_datastores)
            .def("list_kmer_counters",&WorkSpace::list_kmer_counters)
            .def("list_distance_graphs",&WorkSpace::list_distance_graphs)
            .def("dump",&WorkSpace::dump_to_disk,"filename"_a)
            .def("ls",&WorkSpace::ls,"level"_a=0,"recursive"_a=true)
            ;

    py::class_<PerfectMatch>(m,"PerfectMatch", "A perfect match between a read and a node")
            .def(py::init<sgNodeID_t,int32_t,int32_t,uint16_t>(),"node"_a=0,"node_position"_a=0,"read_position"_a=0,"size"_a=0,py::return_value_policy::take_ownership)
            .def_readonly("node",&PerfectMatch::node)
            .def_readonly("node_position",&PerfectMatch::node_position)
            .def_readonly("read_position",&PerfectMatch::read_position)
            .def_readonly("size",&PerfectMatch::size)
            .def("__repr__", [](const PerfectMatch &pm) {
                return "<PerfectMatch: " + std::to_string(pm.read_position) + ":" +
                    std::to_string(pm.read_position+pm.size-1) + " -> " + std::to_string(pm.node) + " @ " +
                    std::to_string(pm.node_position) +":"+ std::to_string(pm.node_position+pm.size-1)+">";
            })
            ;

    py::class_<PerfectMatchesFilter>(m,"PerfectMatchesFilter","Collection of static methods that filter PerfercMatch vectors")
            .def(py::init<WorkSpace&>(),"workspace"_a=0,py::return_value_policy::take_ownership)
            .def("truncate_turnaround",&PerfectMatchesFilter::truncate_turnaround,"matches"_a,py::return_value_policy::take_ownership)
            .def("matches_fw_from_node",&PerfectMatchesFilter::matches_fw_from_node,"node"_a,"matches"_a,py::return_value_policy::take_ownership)
            .def("clean_linear_groups",&PerfectMatchesFilter::clean_linear_groups,"matches"_a,"group_size"_a=5,"small_node_size"_a=500,py::return_value_policy::take_ownership)
            .def("merge_and_sort",&PerfectMatchesFilter::merge_and_sort,"vvmatches"_a,py::return_value_policy::take_ownership)
            ;

    py::class_<PerfectMatchesMergeSorter>(m,"PerfectMatchesMergeSorter","A whole class to merge multiple LRs from a node")
            .def(py::init<WorkSpace&>(),"workspace"_a=0,py::return_value_policy::take_ownership)
            .def("init_from_node",&PerfectMatchesMergeSorter::init_from_node,"node"_a,"lrr"_a,"min_reads"_a=3, "group_size"_a=5,"small_node_size"_a=500)
            .def("find_next_node",&PerfectMatchesMergeSorter::find_next_node,"d"_a=2000,"candidate_percentaje"_a=0.5,"first_percentaje"_a=0.8, "verbose"_a=false)
            .def("advance_reads_to_node",&PerfectMatchesMergeSorter::advance_reads_to_node)
            .def("advance_reads_through_node",&PerfectMatchesMergeSorter::advance_reads_through_node)
            .def("drop_conflictive_reads",&PerfectMatchesMergeSorter::drop_conflictive_reads)
            .def_readwrite("next_node",&PerfectMatchesMergeSorter::next_node)
            .def_readwrite("read_matches",&PerfectMatchesMergeSorter::read_matches)
            .def_readwrite("read_next_match",&PerfectMatchesMergeSorter::read_next_match)
            .def_readwrite("read_dropped_position",&PerfectMatchesMergeSorter::read_dropped_position)
            .def_readwrite("read_last_hit_position",&PerfectMatchesMergeSorter::read_last_hit_position)
            .def_readwrite("out",&PerfectMatchesMergeSorter::out)
            ;

    py::class_<LongReadsRecruiter>(m, "LongReadsRecruiter", "Long reads exact mapper, threader, etd")
            .def(py::init<SequenceDistanceGraph &, const LongReadsDatastore &,uint8_t, uint16_t>(),"sdg"_a,"datastore"_a,"k"_a=25,"f"_a=50,py::return_value_policy::take_ownership)
            .def_readwrite("read_perfect_matches",&LongReadsRecruiter::read_perfect_matches)
            .def_readwrite("node_reads",&LongReadsRecruiter::node_reads)
            .def_readwrite("read_threads",&LongReadsRecruiter::read_threads)
            .def_readwrite("node_threads",&LongReadsRecruiter::node_threads)
            .def_readwrite("read_paths",&LongReadsRecruiter::read_paths)
            .def_readwrite("node_paths",&LongReadsRecruiter::node_paths)
            .def("dump",&LongReadsRecruiter::dump,"filename"_a)
            .def("dump_threads",&LongReadsRecruiter::dump_threads,"filename"_a)
            .def("load",&LongReadsRecruiter::load,"filename"_a)
            .def("load_threads",&LongReadsRecruiter::load_threads,"filename"_a)
            .def("map_old",&LongReadsRecruiter::perfect_mappings,"hit_size"_a=21,"first_read"_a=1,"last_read"_a=0)
            .def("map",&LongReadsRecruiter::map,"hit_size"_a=21,"first_read"_a=1,"last_read"_a=0)
            .def("recruit",&LongReadsRecruiter::recruit_reads,"hit_size"_a=21,"hit_count"_a=1,"first_read"_a=1,"last_read"_a=0)
            .def("recruit_threads",&LongReadsRecruiter::recruit_threads)
            .def("thread_reads",&LongReadsRecruiter::thread_reads,"end_size"_a=500,"end_matches"_a=2)
            .def("simple_thread_reads",&LongReadsRecruiter::simple_thread_reads)
            .def("dg_from_threads",&LongReadsRecruiter::dg_from_threads,"multi_link"_a=false,"remove_duplicated"_a=false,"min_thread_nodes"_a=1)
            .def("rtg_from_threads",&LongReadsRecruiter::rtg_from_threads,"remove_duplicated"_a=false,"min_thread_nodes"_a=1)
            .def("endmatches_to_positions",&LongReadsRecruiter::endmatches_to_positions)
            .def("thread_and_pop",&LongReadsRecruiter::thread_and_pop)
            .def("path_fw",&LongReadsRecruiter::path_fw,"read_id"_a,"node"_a)
            .def("all_paths_fw",&LongReadsRecruiter::all_paths_fw,"node"_a)
            .def("reverse_perfect_matches",&LongReadsRecruiter::reverse_perfect_matches,"matches"_a,"read_size"_a=0)
            ;

    py::class_<NodePosition>(m,"NodePosition", "A node position in a LRR")
            .def(py::init<sgNodeID_t ,int32_t ,int32_t >(),"node"_a,"start"_a,"end"_a,py::return_value_policy::take_ownership)
            .def_readwrite("node",&NodePosition::node)
            .def_readwrite("start",&NodePosition::start)
            .def_readwrite("end",&NodePosition::end)
            .def("__repr__", [](const NodePosition &np) {
                return "<NodePosition: node " + std::to_string(np.node) + " @ "+std::to_string(np.start)+":"+std::to_string(np.end)+">";
            })
            ;


    py::class_<GraphEditor>(m,"GraphEditor", "A graph editor with operation queue")
            .def(py::init<WorkSpace &>(),py::return_value_policy::take_ownership)
            .def("queue_node_deletion",&GraphEditor::queue_node_deletion)
            .def("queue_link_deletion",&GraphEditor::queue_link_deletion)
            .def("queue_node_expansion",&GraphEditor::queue_node_expansion)
            .def("queue_path_detachment",&GraphEditor::queue_path_detachment)
            .def("apply_all",&GraphEditor::apply_all)
            .def("remove_small_components",&GraphEditor::remove_small_components)
            .def("__repr__", [](const GraphEditor &ge) {
                return "<Graph editor: Node deletions queued: " + std::to_string(ge.node_deletion_queue.size()) +
                        ", link deletion queued: " + std::to_string(ge.link_deletion_queue.size()) +
                        ", path detachment queued: " + std::to_string(ge.path_detachment_queue.size()) +
                        ", node expansion queued: " + std::to_string(ge.node_expansion_queue.size()) +
                        ", next_ops: " + std::to_string(ge.next_op);
            })
            ;

    py::class_<LinkageUntangler>(m,"LinkageUntangler", "A linkage untangler")
            .def(py::init<const DistanceGraph &>(),py::return_value_policy::take_ownership)
            .def_readonly("selected_nodes",&LinkageUntangler::selected_nodes)
            .def("select_by_size",&LinkageUntangler::select_by_size, "min_size"_a, "max_size"_a=0)
            .def("select_all",&LinkageUntangler::select_all)
            .def("select_linear_anchors",&LinkageUntangler::select_linear_anchors,"min_links"_a=3,"min_transitive_links"_a=2)
            .def("make_nextselected_linkage",&LinkageUntangler::make_nextselected_linkage,"min_links"_a=3)
            ;

    py::class_<GraphMaker>(m,"GraphMaker","DBG construccion")
            .def(py::init<SequenceDistanceGraph &>(),py::return_value_policy::take_ownership)
            .def("new_graph_from_paired_datastore",&GraphMaker::new_graph_from_paired_datastore)
            .def("new_graph_from_long_datastore",&GraphMaker::new_graph_from_long_datastore)
            ;

    py::class_<GraphContigger>(m,"GraphContigger","Paired end contigger")
            .def(py::init<WorkSpace &>(),py::return_value_policy::take_ownership)
            .def("reconnect_tips",&GraphContigger::reconnect_tips,"datastore"_a,"min_links"_a=3)
            .def("clip_tips",&GraphContigger::clip_tips,"size"_a,"rounds"_a=10)
            .def("pop_bubbles",&GraphContigger::pop_bubbles,"datastore"_a,"bubble_size"_a=200,"min_support"_a=4,"max_noise"_a=3,"snr"_a=10)
            .def("solve_canonical_repeats_with_single_paths",&GraphContigger::solve_canonical_repeats_with_single_paths,"datastore"_a,"min_support"_a,"max_noise"_a,"snr"_a=10, "join_unitigs"_a=true, "dry_run"_a=false, "verbose"_a=false)
            .def("solve_canonical_repeats_with_paired_paths",&GraphContigger::solve_canonical_repeats_with_paired_paths,"datastore"_a,"min_support"_a,"max_noise"_a,"snr"_a=10, "join_unitigs"_a=true, "dry_run"_a=false, "verbose"_a=false)
            .def("solve_canonical_repeats_with_long_reads",&GraphContigger::solve_canonical_repeats_with_long_reads,"recruiter"_a,"max_kci"_a,"min_support"_a,"max_noise"_a,"snr"_a=10)
            .def("solve_bubble",&GraphContigger::solve_bubble, "t"_a, "s"_a, "ge"_a)
            .def("solve_repeat",&GraphContigger::solve_repeat, "t"_a, "s"_a, "ge"_a)
            .def("solve_tip",&GraphContigger::solve_tip, "t"_a, "s"_a, "ge"_a)
            .def("solve_unclassified",&GraphContigger::solve_unclassified, "t"_a, "s"_a, "ge"_a)
            .def("solve_canonical_repeat",&GraphContigger::solve_canonical_repeat, "ge"_a, "nv"_a, "peds"_a, "min_support"_a=5, "max_noise"_a=10, "snr"_a=10, "verbose"_a=false)
            .def("clip_tip",&GraphContigger::clip_tip, "ge"_a, "nv"_a, "peds"_a, "min_support"_a=5, "max_noise"_a=10, "snr"_a=10, "verbose"_a=false)
            .def("pop_error_bubble",&GraphContigger::pop_error_bubbble, "ge"_a, "nv1"_a, "nv2"_a, "peds"_a, "min_support"_a=5, "max_noise"_a=10, "snr"_a=10, "verbose"_a=false)
            .def("solve_all_tangles",&GraphContigger::solve_all_tangles, "ws"_a, "peds"_a, "fsize"_a=220, "fminkci"_a=-1, "fmaxkci"_a=-1, "apply"_a=false)
            .def("solve_all_canonical",&GraphContigger::solve_all_canonical, "ge"_a, "peds"_a, "size"_a=1000, "apply"_a=false)
            .def("clip_all_tips",&GraphContigger::clip_all_tips, "ge"_a, "peds"_a, "size"_a=1000, "apply"_a=false)
            .def("pop_all_error_bubbles",&GraphContigger::pop_all_error_bubbles, "ge"_a, "peds"_a, "size"_a=1000, "apply"_a=false)
            .def("contig_reduction_to_unique_kmers",&GraphContigger::contig_reduction_to_unique_kmers, "ws"_a, "kmer_counter"_a, "kmer_count"_a, "min_cov"_a, "max_cov"_a, "max_run_size"_a=100)
            .def("group_nodes",&GraphContigger::group_nodes, "peds"_a)
            ;

    py::class_<LinkageMaker>(m, "LinkageMaker", "Makes linkage")
            .def(py::init<SequenceDistanceGraph &>())
            .def("select_by_size",&LinkageMaker::select_by_size,"min_size"_a,"max_size"_a=0)
            .def("make_longreads_multilinkage",&LinkageMaker::make_longreads_multilinkage,"lorm"_a, "min_map_size"_a=1000, "min_map_id"_a=.1, "real_read_size"_a=true, "unmapped_end"_a=1000)
            ;

    py::class_<PathFinder>(m,"PathFinder","PathFinder")
            .def(py::init<WorkSpace &,sgNodeID_t, sgNodeID_t, uint8_t>(),"ws"_a,"n1"_a,"n2"_a,"k"_a,py::return_value_policy::take_ownership)
            .def("index_paths",&PathFinder::index_paths,"paths"_a)
            .def("seq_to_pathpos",&PathFinder::seq_to_pathpos,"path_id"_a,"seq"_a)
            .def("load_lrseqs",&PathFinder::load_lrseqs,"dg"_a,"lrr"_a,"ovl_extension"_a=300)
            .def("lrseqs_as_fasta",&PathFinder::lrseqs_as_fasta)
            .def("index_seqs",&PathFinder::index_seqs)
            ;

    py::class_<PFScoredPath>(m,"PFScoredPath","PFScoredPath")
            .def(py::init<PathFinder &, sgNodeID_t ,sgNodeID_t>(),"pf"_a,"from"_a,"to"_a,py::return_value_policy::take_ownership)
            .def("find_hits",&PFScoredPath::find_hits)
            .def_readonly("read_hitpos",&PFScoredPath::read_hitpos,py::return_value_policy::reference)
            .def_readwrite("path",&PFScoredPath::path,py::return_value_policy::reference)
            .def("score",&PFScoredPath::score)
            ;

    py::class_<Strider>(m,"Strider","Strider")
            .def(py::init<WorkSpace &>(),"ws"_a,py::return_value_policy::take_ownership)
            .def_readonly_static("logo",&Strider::logo)
            .def_readwrite("experimental_striding",&Strider::experimental_striding)
            .def_readonly("routes_fw",&Strider::routes_fw)
            .def_readonly("routes_bw",&Strider::routes_bw)
            .def_readonly("links_fw",&Strider::links_fw)
            .def_readonly("links_bw",&Strider::links_bw)
            .def_readonly("is_anchor",&Strider::is_anchor)
            .def("dump",&Strider::dump,"filename"_a)
            .def("load",&Strider::load,"filename"_a)
            .def("route_vs_readpaths_stats",&Strider::route_vs_readpaths_stats)
            .def("join_stride_single_strict_from_all",&Strider::join_stride_single_strict_from_all)
            .def("stride_from_anchors",&Strider::stride_from_anchors,"min_size"_a=1,"min_kci"_a=.5,"max_kci"_a=1.5)
            .def("link_from_anchors",&Strider::link_from_anchors,"min_size"_a=1,"min_kci"_a=.5,"max_kci"_a=1.5,"d"_a=2000, "min_reads"_a=3, "group_size"_a=5, "small_node_size"_a=500,"candidate_percentaje"_a=0.5,"first_percentaje"_a=0.8)
            .def("add_datastore",py::overload_cast<const PairedReadsDatastore &>(&Strider::add_datastore))
            .def("add_datastore",py::overload_cast<const LongReadsRecruiter &>(&Strider::add_datastore))
            .def("stride_out",&Strider::stride_out,"node"_a,py::return_value_policy::take_ownership)
            .def("stride_single_strict",&Strider::stride_single_strict,"node"_a,"min_reads"_a=3,"max_noise"_a=.1,py::return_value_policy::take_ownership)
            .def("stride_out_in_order",&Strider::stride_out_in_order,"node"_a,"use_pair"_a=true,"collapse_pair"_a=true,"verbose"_a=false,py::return_value_policy::take_ownership)
            .def("link_out_by_lr",&Strider::link_out_by_lr,"node"_a,"d"_a=2000, "min_reads"_a=3, "group_size"_a=5, "small_node_size"_a=500, "candidate_percentaje"_a=0.5, "first_percentaje"_a=0.8, "verbose"_a=false, py::return_value_policy::take_ownership)
            ;

    py::class_<GraphPatcher>(m,"GraphPatcher","GraphPatcher")
            .def(py::init<WorkSpace &>(),"ws"_a,py::return_value_policy::take_ownership)
            .def("find_tips_to_reconnect",&GraphPatcher::find_tips_to_reconnect)
            .def("collapse_reconnection_groups",&GraphPatcher::collapse_reconnection_groups)
            .def("patch_reconnection_groups",&GraphPatcher::patch_reconnection_groups)
            .def("create_patch",&GraphPatcher::create_patch,"reconnection_group"_a)
            .def_readonly("reconnection_groups",&GraphPatcher::reconnection_groups);

    py::class_<LineFiller>(m, "LineFiller", "LineFiller")
            .def(py::init<WorkSpace &, LongReadsRecruiter&>(), "ws"_a, "lrr"_a)
            .def("line_fill", &LineFiller::line_fill, "anchor_path"_a)
            .def("fill_all_paths", &LineFiller::fill_all_paths, "lines"_a)
            .def_readonly("final_lines", &LineFiller::final_lines);


    py::class_<CountFilter>(m,"CountFilter","CountFilter")
            .def(py::init<std::string, std::string , int , std::vector<std::string>, std::vector<int>>(),"kcname"_a, "filter_count_name"_a, "filter_count_max"_a, "value_count_names"_a, "value_count_mins"_a)
            .def("get_pattern",&CountFilter::get_pattern,"nv"_a);

    py::class_<LocalOrder>(m,"LocalOrder","LocalOrder")
            .def(py::init<>())
            .def(py::init<const std::vector<sgNodeID_t> &>(),"nodes"_a)
            .def("get_node_position",&LocalOrder::get_node_position,"node_id"_a)
            .def("place_node",&LocalOrder::place_node,"rtg"_a,"nid"_a, "node_coordinates"_a=std::unordered_map<sgNodeID_t,int64_t>(), "max_hops"_a=3, "min_links"_a=3)
            .def("add_placed_nodes",&LocalOrder::add_placed_nodes,"placed_nodes"_a,"update_current"_a=false)
            .def("remove_node",&LocalOrder::remove_node,"nid"_a)
            .def("as_signed_nodes",&LocalOrder::as_signed_nodes)
            .def("as_thread",&LocalOrder::as_thread,"dg"_a)
            .def("size",&LocalOrder::size)
            .def("reverse",&LocalOrder::reverse)
            .def("thread_order_crosses",&LocalOrder::thread_order_crosses,"thread")
            .def("merge",&LocalOrder::merge, "other"_a, "max_overhang"_a=4, "min_shared_perc"_a=.5, "min_shared"_a=20, "max_disordered_perc"_a=.02)
            .def_readwrite("node_positions",&LocalOrder::node_positions)
            .def_readwrite("node_coordinates",&LocalOrder::node_coordinates);

    py::class_<RTGCluster>(m, "RTGCluster", "RTGCluster")
            .def(py::init<const ReadThreadsGraph &, int, int, int, float>(),"","rtg"_a,"p"_a,"q"_a,"min_node_threads"_a,"min_node_happiness"_a, py::return_value_policy::take_ownership)
            .def(py::init<const RTGCluster & >(),"","rtgcluster"_a, py::return_value_policy::take_ownership)
            .def("is_node_happy",&RTGCluster::is_node_happy,"nid"_a)
            .def("new_happy_nodes",&RTGCluster::new_happy_nodes)
            .def("add_node",&RTGCluster::add_node,"nid"_a)
            .def("is_thread_happy",&RTGCluster::is_thread_happy,"tid"_a)
            .def("new_happy_threads",&RTGCluster::new_happy_threads)
            .def("add_thread",&RTGCluster::add_thread,"tid"_a)
            .def("grow",&RTGCluster::grow,"steps"_a=UINT64_MAX)
            .def_readonly("nodes",&RTGCluster::nodes)
            .def_readonly("threads",&RTGCluster::threads)
            .def_readonly("node_total_threads",&RTGCluster::node_total_threads)
            .def_readonly("node_happiness_threads",&RTGCluster::node_happiness_threads)
            .def_readonly("node_threads",&RTGCluster::node_threads)
            .def_readonly("thread_nodes",&RTGCluster::thread_nodes)
            ;

    py::class_<RTGClassifier>(m, "RTGClassifier", "RTGClassifier")
            .def(py::init<const ReadThreadsGraph &, int, float, int, int>(),"","rtg"_a,"min_node_threads"_a=2,"min_node_percentage"_a=.6,"thread_p"_a=10,"thread_q"_a=12, py::return_value_policy::take_ownership)
            .def("compute_node_class",&RTGClassifier::compute_node_class,"nid"_a)
            .def("switch_node_class",&RTGClassifier::switch_node_class,"nid"_a,"c"_a)
            .def("compute_thread_class",&RTGClassifier::compute_thread_class,"tid"_a)
            .def("switch_thread_class",&RTGClassifier::switch_thread_class,"tid"_a,"c"_a)
            .def("propagate",&RTGClassifier::propagate,"steps"_a=UINT64_MAX,"verbose"_a=false)
            .def_readonly("node_threads",&RTGClassifier::node_threads)
            .def_readonly("thread_nodes",&RTGClassifier::thread_nodes)
            .def_readonly("node_class",&RTGClassifier::node_class)
            .def_readonly("thread_class",&RTGClassifier::thread_class)
            .def_readonly("nodes_to_evaluate",&RTGClassifier::nodes_to_evaluate)
            .def_readonly("threads_to_evaluate",&RTGClassifier::threads_to_evaluate)
            ;

    m.def("str_to_kmers",&sdglib::str_to_kmers,py::return_value_policy::take_ownership);
    m.def("str_rc",&sdglib::str_rc,py::return_value_policy::take_ownership);
    m.def("count_kmers_as_graph_nodes",&BatchKmersCounter::countKmersToGraphNodes,py::return_value_policy::take_ownership,"sdg"_a,"peds"_a,"k"_a,"min_coverage"_a, "max_coverage"_a, "num_batches"_a);
}
