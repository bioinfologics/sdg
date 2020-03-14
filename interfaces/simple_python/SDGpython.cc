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

namespace py = pybind11;
using namespace py::literals;


PYBIND11_MAKE_OPAQUE(std::vector<std::vector<PerfectMatch>>);
PYBIND11_MAKE_OPAQUE(std::vector<std::vector<uint64_t>>);
PYBIND11_MAKE_OPAQUE(std::vector<std::vector<NodePosition>>);

PYBIND11_MODULE(SDGpython, m) {

    py::bind_vector<std::vector<uint64_t>>(m, "VectorU64");
    py::bind_vector<std::vector<std::vector<uint64_t>>>(m, "VectorVectorU64");
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
            .export_values();

    py::enum_<KmerCountMode>(m,"KmerCountMode")
            .value("Canonical",KmerCountMode::Canonical)
            .value("NonCanonical",KmerCountMode::NonCanonical)
            .export_values();

    py::class_<Support>(m,"Support","Support for an operation or piece of data")
            .def(py::init<SupportType,uint16_t,uint64_t>(),"","support"_a=SupportType::Undefined,"index"_a=0,"id"_a=0)
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
            .def(py::init<const SequenceDistanceGraph &>(),"","sdg"_a)
            .def("add_link",&DistanceGraph::add_link,"source"_a,"dest"_a,"distance"_a,"support"_a=Support() )
            .def("disconnect_node",&DistanceGraph::disconnect_node)
            .def("remove_link",py::overload_cast<sgNodeID_t , sgNodeID_t >(&DistanceGraph::remove_link))
            .def("remove_link",py::overload_cast<sgNodeID_t , sgNodeID_t, int32_t, Support >(&DistanceGraph::remove_link))
            .def("fw_neighbours_by_distance",&DistanceGraph::fw_neighbours_by_distance,"node_id"_a,"min_links"_a=3)
            .def("get_nodeview",&DistanceGraph::get_nodeview)
            .def("get_all_nodeviews",&DistanceGraph::get_all_nodeviews,"include_disconnected"_a=true,"Returns a vector with NodeViews for all active nodes",py::return_value_policy::move)
            .def_readwrite("name",&DistanceGraph::name)
            .def("write_to_gfa1",&DistanceGraph::write_to_gfa1,"filename"_a,"selected_nodes"_a=std::vector<sgNodeID_t>(),"depths"_a=std::vector<sgNodeID_t>())
            .def("write_to_gfa2",&DistanceGraph::write_to_gfa2)
            ;

    py::class_<SequenceDistanceGraph,DistanceGraph>(m, "SequenceDistanceGraph", "A Sequence Distance Graph")
            .def("get_node_size",&SequenceDistanceGraph::get_node_size)
            .def("add_node",py::overload_cast<std::string>(&SequenceDistanceGraph::add_node))
            .def("remove_node",&SequenceDistanceGraph::remove_node)
            .def("join_all_unitigs",&SequenceDistanceGraph::join_all_unitigs)
            .def("load_from_gfa",&SequenceDistanceGraph::load_from_gfa)
            .def("load_from_fasta",&SequenceDistanceGraph::load_from_fasta)
            ;

    py::class_<PairedReadsDatastore>(m, "PairedReadsDatastore", "A Paired Reads Datastore")
            .def("size",&PairedReadsDatastore::size)
            .def("get_read_sequence",&PairedReadsDatastore::get_read_sequence)
            .def_readwrite("mapper",&PairedReadsDatastore::mapper)
            ;

    py::class_<PairedReadsMapper>(m, "PairedReadsMapper", "A Paired Reads Mapper")
            .def("path_reads63",&PairedReadsMapper::path_reads63)
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
            .def("count_spectra",&KmerCounter::count_spectra,"name"_a,"max_freq"_a=1000,"unique_in_graph"_a=true)
            .def("update_graph_counts",&KmerCounter::update_graph_counts)
            ;

    py::class_<WorkSpace>(m, "WorkSpace", "A full SDG WorkSpace")
            .def(py::init<>(),"")
            .def(py::init<const std::string &>(),"Load from disk","filename"_a)
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
            .def_readonly("node",&PerfectMatch::node)
            .def_readonly("node_position",&PerfectMatch::node_position)
            .def_readonly("read_position",&PerfectMatch::read_position)
            .def_readonly("size",&PerfectMatch::size)
            ;
    py::bind_vector<std::vector<PerfectMatch>>(m, "VectorPerfectMatch");
    py::bind_vector<std::vector<std::vector<PerfectMatch>>>(m, "VectorVectorPerfectMatch");

    py::class_<LongReadsRecruiter>(m, "LongReadsRecruiter", "A full SDG WorkSpace")
            .def(py::init<const SequenceDistanceGraph &, const LongReadsDatastore &,uint8_t, uint16_t>(),"sdg"_a,"datastore"_a,"k"_a=25,"f"_a=50)
            .def_readwrite("read_perfect_matches",&LongReadsRecruiter::read_perfect_matches)
            .def_readwrite("node_reads",&LongReadsRecruiter::node_reads)
            .def_readwrite("read_threads",&LongReadsRecruiter::read_threads)
            .def("dump",&LongReadsRecruiter::dump,"filename"_a)
            .def("load",&LongReadsRecruiter::load,"filename"_a)
            .def("map",&LongReadsRecruiter::perfect_mappings,"hit_size"_a=21,"first_read"_a=1,"last_read"_a=0)
            .def("recruit",&LongReadsRecruiter::recruit_reads,"hit_size"_a=21,"hit_count"_a=1,"first_read"_a=1,"last_read"_a=0)
            .def("thread_reads",&LongReadsRecruiter::thread_reads,"end_size"_a=500,"end_matches"_a=2)
            .def("dg_from_threads",&LongReadsRecruiter::dg_from_threads,"multi_link"_a=false)
            .def("endmatches_to_positions",&LongReadsRecruiter::endmatches_to_positions)
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
    py::bind_vector<std::vector<std::vector<NodePosition>>>(m, "VectorVectorNodePosition");

    py::class_<GraphEditor>(m,"GraphEditor", "A graph editor with operation queue")
            .def(py::init<WorkSpace &>())
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
            ;

    py::class_<GraphContigger>(m,"GraphContigger","Paired end contigger")
            .def(py::init<WorkSpace &>())
            .def("clip_tips",&GraphContigger::clip_tips,"size"_a,"rounds"_a=10)
            .def("pop_bubbles",&GraphContigger::pop_bubbles,"datastore"_a,"bubble_size"_a=200,"min_support"_a=4,"max_noise"_a=3,"snr"_a=10)
            .def("solve_canonical_repeats_with_single_paths",&GraphContigger::solve_canonical_repeats_with_single_paths,"datastore"_a,"min_support"_a,"max_noise"_a,"snr"_a=10)
            .def("solve_canonical_repeats_with_paired_paths",&GraphContigger::solve_canonical_repeats_with_paired_paths,"datastore"_a,"min_support"_a,"max_noise"_a,"snr"_a=10)
            ;
}