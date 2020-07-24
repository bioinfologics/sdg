//
// Created by Bernardo Clavijo (EI) on 13/02/2020.
//

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>
#include <sdglib/graph/SequenceDistanceGraph.hpp>
#include <sdglib/workspace/WorkSpace.hpp>
#include <sdglib/views/NodeView.hpp>
#include <sdglib/mappers/LongReadsRecruiter.hpp>
#include <sdglib/processors/GraphEditor.hpp>
#include <sdglib/processors/LinkageUntangler.hpp>
#include <sdglib/processors/GraphContigger.hpp>
#include <sdglib/processors/GraphMaker.hpp>
#include <sdglib/processors/GraphPatcher.hpp>
#include <sdglib/processors/PathFinder.hpp>
#include <sdglib/processors/Strider.hpp>
#include <sdglib/processors/LineFiller.h>

namespace py = pybind11;
using namespace py::literals;


PYBIND11_MAKE_OPAQUE(std::vector<std::vector<PerfectMatch>>);
PYBIND11_MAKE_OPAQUE(std::vector<std::vector<uint64_t>>);
PYBIND11_MAKE_OPAQUE(std::vector<std::vector<int64_t>>);
PYBIND11_MAKE_OPAQUE(std::vector<std::vector<NodePosition>>);
PYBIND11_MAKE_OPAQUE(std::vector<ReadPath>);
PYBIND11_MAKE_OPAQUE(std::vector<bool>);
PYBIND11_MAKE_OPAQUE(std::vector<std::vector<Link>>);

PYBIND11_MODULE(SDGpython, m) {

    py::bind_vector<std::vector<uint64_t>>(m, "VectorU64");
    py::bind_vector<std::vector<std::vector<uint64_t>>>(m, "VectorVectorU64");
    py::bind_vector<std::vector<int64_t>>(m, "Vector64");
    py::bind_vector<std::vector<std::vector<int64_t>>>(m, "VectorVector64");
    py::bind_vector<std::vector<ReadPath>>(m, "VectorReadPath");
    py::bind_vector<std::vector<PerfectMatch>>(m, "VectorPerfectMatch");
    py::bind_vector<std::vector<std::vector<PerfectMatch>>>(m, "VectorVectorPerfectMatch");
    py::bind_vector<std::vector<std::vector<Link>>>(m, "VectorVectorLink");
    py::bind_vector<std::vector<std::vector<NodePosition>>>(m, "VectorVectorNodePosition");
    py::bind_vector<std::vector<bool>>(m, "VectorBool");

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
            .def(py::init<SupportType,uint16_t,int64_t>(),"","support"_a=SupportType::Undefined,"index"_a=0,"id"_a=0)
            .def_readonly("type",&Support::type)
            .def_readonly("index",&Support::index)
            .def_readonly("id",&Support::id)
            ;
    py::class_<Node>(m, "Node", "A node in a Sequence Distance Graph")
            .def(py::init<const std::string &>(),"","sequence"_a)
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
            .def(py::init<sgNodeID_t, sgNodeID_t, int32_t, Support>(),"","src"_a,"dest"_a,"dist"_a,"support"_a=Support())
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
            .def("sequence",&NodeView::sequence)
            .def("size",&NodeView::size)
            .def("prev",&NodeView::prev,"A list with LinkViews to previous nodes",py::return_value_policy::move)
            .def("next",&NodeView::next,"A list with LinkViews to next nodes",py::return_value_policy::move)
            .def("parallels",&NodeView::parallels,"A list with NodeViews of parallel nodes",py::return_value_policy::move)
            .def("kmer_coverage",py::overload_cast<std::string, std::string>(&NodeView::kmer_coverage, py::const_))
            .def("kci",&NodeView::kci)
            .def("__eq__", &NodeView::operator==, py::is_operator())
            .def("__repr__",
                 [](const NodeView &nv) {
                     return "<NodeView: node " + std::to_string(nv.node_id()) + " in graph " + nv.graph().name + ">";
                 })
            ;

    py::class_<DistanceGraph>(m, "DistanceGraph", "A Distance Graph")
            .def_property_readonly("sdg", [] (const DistanceGraph &dg) {return &dg.sdg;},py::return_value_policy::reference)
            .def(py::init<SequenceDistanceGraph &, const std::string &>(),"","sdg"_a,"name"_a="unnamed")
            .def("add_link",&DistanceGraph::add_link,"source"_a,"dest"_a,"distance"_a,"support"_a=Support() )
            .def("disconnect_node",&DistanceGraph::disconnect_node)
            .def("remove_link",py::overload_cast<sgNodeID_t , sgNodeID_t >(&DistanceGraph::remove_link))
            .def("remove_link",py::overload_cast<sgNodeID_t , sgNodeID_t, int32_t, Support >(&DistanceGraph::remove_link))
            .def("remove_transitive_links",&DistanceGraph::remove_transitive_links,"radius"_a=10)
            .def("fw_neighbours_by_distance",&DistanceGraph::fw_neighbours_by_distance,"node_id"_a,"min_links"_a=3)
            .def("fw_reached_nodes",&DistanceGraph::fw_reached_nodes,"node_id"_a,"radius"_a=10)
            .def("get_nodeview",&DistanceGraph::get_nodeview)
            .def("get_all_nodeviews",&DistanceGraph::get_all_nodeviews,"both_directions"_a=false,"include_disconnected"_a=true,"Returns a vector with NodeViews for all active nodes",py::return_value_policy::move)
            .def_readwrite("name",&DistanceGraph::name)
            .def("write_to_gfa1",&DistanceGraph::write_to_gfa1,"filename"_a,"selected_nodes"_a=std::vector<sgNodeID_t>(),"depths"_a=std::vector<sgNodeID_t>())
            .def("write_to_gfa2",&DistanceGraph::write_to_gfa2)
            .def("stats_by_kci",&DistanceGraph::stats_by_kci)
            .def("are_connected",&DistanceGraph::are_connected,"node1"_a, "node2"_a)
            .def("find_all_paths_between",&DistanceGraph::find_all_paths_between,"from"_a,"to"_a,"max_size"_a,"max_nodes"_a,"abort_on_loops"_a)
            .def("get_all_lines",&DistanceGraph::get_all_lines,"min_nodes"_a, "min_total_size"_a=0)
            ;

    py::class_<SequenceDistanceGraph,DistanceGraph>(m, "SequenceDistanceGraph", "A Sequence Distance Graph")
            .def("get_node_size",&SequenceDistanceGraph::get_node_size)
            .def("add_node",py::overload_cast<std::string>(&SequenceDistanceGraph::add_node))
            .def("remove_node",&SequenceDistanceGraph::remove_node)
            .def("join_all_unitigs",&SequenceDistanceGraph::join_all_unitigs)
            .def("load_from_gfa",&SequenceDistanceGraph::load_from_gfa)
            .def("load_from_bcalm",&SequenceDistanceGraph::load_from_bcalm,"filename"_a,"k"_a)
            .def("load_from_fasta",&SequenceDistanceGraph::load_from_fasta)
            ;

    py::class_<SequenceDistanceGraphPath>(m, "SequenceDistanceGraphPath", "SequenceDistanceGraphPath")
            .def(py::init<const SequenceDistanceGraph &,const std::vector<sgNodeID_t> >(),"","sdg"_a,"nodes"_a)
            .def("sequence", &SequenceDistanceGraphPath::sequence)
            .def_readwrite("nodes", &SequenceDistanceGraphPath::nodes)
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
            .def("path_reads63",&PairedReadsMapper::path_reads63)
            .def("path_reads",&PairedReadsMapper::path_reads)
            .def_readonly("paths_in_node",&PairedReadsMapper::paths_in_node)
            .def_readonly("read_paths",&PairedReadsMapper::read_paths)
            .def("path_fw",&PairedReadsMapper::path_fw,"read_id"_a,"node"_a,"use_pair"_a=true,"collapse_pair"_a=true)
            .def("all_paths_fw",&PairedReadsMapper::all_paths_fw,"node"_a,"use_pair"_a=true,"collapse_pair"_a=true)
            .def("dump_readpaths",&PairedReadsMapper::dump_readpaths)
            .def("load_readpaths",&PairedReadsMapper::load_readpaths)
            ;
    py::class_<LinkedReadsDatastore>(m, "LinkedReadsDatastore", "A Linked Reads Datastore")
            .def("size",&LinkedReadsDatastore::size)
            .def("get_read_sequence",&LinkedReadsDatastore::get_read_sequence)
            ;
    py::class_<LongReadsDatastore>(m, "LongReadsDatastore", "A Long Reads Datastore")
            .def("size",&LongReadsDatastore::size)
            .def("get_read_sequence",&LongReadsDatastore::get_read_sequence)
            .def("get_read_size",[](const LongReadsDatastore &ds,__uint64_t rid) {
                return ds.read_to_fileRecord[rid].record_size;
            })
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

    py::class_<LongReadsRecruiter>(m, "LongReadsRecruiter", "A full SDG WorkSpace")
            .def(py::init<SequenceDistanceGraph &, const LongReadsDatastore &,uint8_t, uint16_t>(),"sdg"_a,"datastore"_a,"k"_a=25,"f"_a=50)
            .def_readwrite("read_perfect_matches",&LongReadsRecruiter::read_perfect_matches)
            .def_readwrite("node_reads",&LongReadsRecruiter::node_reads)
            .def_readwrite("read_threads",&LongReadsRecruiter::read_threads)
            .def_readwrite("node_threads",&LongReadsRecruiter::node_threads)
            .def_readwrite("read_paths",&LongReadsRecruiter::read_paths)
            .def_readwrite("node_paths",&LongReadsRecruiter::node_paths)
            .def("dump",&LongReadsRecruiter::dump,"filename"_a)
            .def("load",&LongReadsRecruiter::load,"filename"_a)
            .def("map_old",&LongReadsRecruiter::perfect_mappings,"hit_size"_a=21,"first_read"_a=1,"last_read"_a=0)
            .def("map",&LongReadsRecruiter::map,"hit_size"_a=21,"first_read"_a=1,"last_read"_a=0)
            .def("recruit",&LongReadsRecruiter::recruit_reads,"hit_size"_a=21,"hit_count"_a=1,"first_read"_a=1,"last_read"_a=0)
            .def("recruit_threads",&LongReadsRecruiter::recruit_threads)
            .def("thread_reads",&LongReadsRecruiter::thread_reads,"end_size"_a=500,"end_matches"_a=2)
            .def("dg_from_threads",&LongReadsRecruiter::dg_from_threads,"multi_link"_a=false)
            .def("endmatches_to_positions",&LongReadsRecruiter::endmatches_to_positions)
            .def("thread_and_pop",&LongReadsRecruiter::thread_and_pop)
            .def("path_fw",&LongReadsRecruiter::path_fw,"read_id"_a,"node"_a)
            .def("all_paths_fw",&LongReadsRecruiter::all_paths_fw,"node"_a)
            .def("reverse_perfect_matches",&LongReadsRecruiter::reverse_perfect_matches,"matches"_a,"read_size"_a=0)
            ;

    py::class_<NodePosition>(m,"NodePosition", "A node position in a LRR")
            .def(py::init<sgNodeID_t ,int32_t ,int32_t >(),"node"_a,"start"_a,"end"_a)
            .def_readwrite("node",&NodePosition::node)
            .def_readwrite("start",&NodePosition::start)
            .def_readwrite("end",&NodePosition::end)
            .def("__repr__", [](const NodePosition &np) {
                return "<NodePosition: node " + std::to_string(np.node) + " @ "+std::to_string(np.start)+":"+std::to_string(np.end)+">";
            })
            ;


    py::class_<GraphEditor>(m,"GraphEditor", "A graph editor with operation queue")
            .def(py::init<WorkSpace &>())
            .def("queue_node_deletion",&GraphEditor::queue_node_deletion)
            .def("queue_link_deletion",&GraphEditor::queue_link_deletion)
            .def("queue_node_expansion",&GraphEditor::queue_node_expansion)
            .def("queue_path_detachment",&GraphEditor::queue_path_detachment)
            .def("apply_all",&GraphEditor::apply_all)
            .def("remove_small_components",&GraphEditor::remove_small_components)
            ;

    py::class_<LinkageUntangler>(m,"LinkageUntangler", "A linkage untangler")
            .def(py::init<const DistanceGraph &>())
            .def_readonly("selected_nodes",&LinkageUntangler::selected_nodes)
            .def("select_all",&LinkageUntangler::select_all)
            .def("select_linear_anchors",&LinkageUntangler::select_linear_anchors,"min_links"_a=3,"min_transitive_links"_a=2)
            .def("make_nextselected_linkage",&LinkageUntangler::make_nextselected_linkage,"min_links"_a=3)
            ;

    py::class_<GraphMaker>(m,"GraphMaker","DBG construccion")
            .def(py::init<SequenceDistanceGraph &>())
            .def("new_graph_from_paired_datastore",&GraphMaker::new_graph_from_paired_datastore)
            .def("new_graph_from_long_datastore",&GraphMaker::new_graph_from_long_datastore)
            ;

    py::class_<GraphContigger>(m,"GraphContigger","Paired end contigger")
            .def(py::init<WorkSpace &>())
            .def("reconnect_tips",&GraphContigger::reconnect_tips,"datastore"_a,"min_links"_a=3)
            .def("clip_tips",&GraphContigger::clip_tips,"size"_a,"rounds"_a=10)
            .def("pop_bubbles",&GraphContigger::pop_bubbles,"datastore"_a,"bubble_size"_a=200,"min_support"_a=4,"max_noise"_a=3,"snr"_a=10)
            .def("solve_canonical_repeats_with_single_paths",&GraphContigger::solve_canonical_repeats_with_single_paths,"datastore"_a,"min_support"_a,"max_noise"_a,"snr"_a=10, "join_unitigs"_a=true, "dry_run"_a=false, "verbose"_a=false)
            .def("solve_canonical_repeats_with_paired_paths",&GraphContigger::solve_canonical_repeats_with_paired_paths,"datastore"_a,"min_support"_a,"max_noise"_a,"snr"_a=10, "join_unitigs"_a=true, "dry_run"_a=false, "verbose"_a=false)
            .def("solve_canonical_repeats_with_long_reads",&GraphContigger::solve_canonical_repeats_with_long_reads,"recruiter"_a,"max_kci"_a,"min_support"_a,"max_noise"_a,"snr"_a=10)
            ;

    py::class_<PathFinder>(m,"PathFinder","PathFinder")
            .def(py::init<WorkSpace &,sgNodeID_t, sgNodeID_t, uint8_t>(),"ws"_a,"n1"_a,"n2"_a,"k"_a)
            .def("index_paths",&PathFinder::index_paths,"paths"_a)
            .def("seq_to_pathpos",&PathFinder::seq_to_pathpos,"path_id"_a,"seq"_a)
            .def("load_lrseqs",&PathFinder::load_lrseqs,"dg"_a,"lrr"_a,"ovl_extension"_a=300)
            .def("lrseqs_as_fasta",&PathFinder::lrseqs_as_fasta)
            .def("index_seqs",&PathFinder::index_seqs)
            ;

    py::class_<PFScoredPath>(m,"PFScoredPath","PFScoredPath")
            .def(py::init<PathFinder &, sgNodeID_t ,sgNodeID_t>(),"pf"_a,"from"_a,"to"_a)
            .def("find_hits",&PFScoredPath::find_hits)
            .def_readonly("read_hitpos",&PFScoredPath::read_hitpos,py::return_value_policy::reference)
            .def_readwrite("path",&PFScoredPath::path,py::return_value_policy::reference)
            .def("score",&PFScoredPath::score)
            ;

    py::class_<Strider>(m,"Strider","Strider")
            .def(py::init<WorkSpace &>(),"ws"_a)
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
            .def("stride_from_anchors",&Strider::stride_from_anchors,"min_size"_a=1,"min_kci"_a=.5,"max_kci"_a=1.5)
            .def("link_from_anchors",&Strider::link_from_anchors,"min_size"_a=1,"min_kci"_a=.5,"max_kci"_a=1.5,"d"_a=2000, "min_reads"_a=3, "group_size"_a=5, "small_node_size"_a=500,"candidate_percentaje"_a=0.5,"first_percentaje"_a=0.8)
            .def("add_datastore",py::overload_cast<const PairedReadsDatastore &>(&Strider::add_datastore))
            .def("add_datastore",py::overload_cast<const LongReadsRecruiter &>(&Strider::add_datastore))
            .def("stride_out",&Strider::stride_out,"node"_a,py::return_value_policy::take_ownership)
            .def("stride_out_in_order",&Strider::stride_out_in_order,"node"_a,"use_pair"_a=true,"collapse_pair"_a=true,"verbose"_a=false,py::return_value_policy::take_ownership)
            .def("link_out_by_lr",&Strider::link_out_by_lr,"node"_a,"d"_a=2000, "min_reads"_a=3, "group_size"_a=5, "small_node_size"_a=500, "candidate_percentaje"_a=0.5, "first_percentaje"_a=0.8, "verbose"_a=false, py::return_value_policy::take_ownership);

    py::class_<GraphPatcher>(m,"GraphPatcher","GraphPatcher")
            .def(py::init<WorkSpace &>(),"ws"_a)
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

    m.def("str_to_kmers",&sdglib::str_to_kmers,py::return_value_policy::take_ownership);
    m.def("str_rc",&sdglib::str_rc,py::return_value_policy::take_ownership);
}
