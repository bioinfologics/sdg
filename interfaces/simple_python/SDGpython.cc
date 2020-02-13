//
// Created by Bernardo Clavijo (EI) on 13/02/2020.
//

#include <deps/pybind11/pybind11.h>
#include <deps/pybind11/stl.h>
#include <deps/pybind11/stl_bind.h>
#include <sdglib/graph/SequenceDistanceGraph.hpp>
#include <sdglib/workspace/WorkSpace.hpp>
#include <sdglib/views/NodeView.hpp>
#include <sdglib/mappers/LongReadsRecruiter.hpp>

namespace py = pybind11;
using namespace py::literals;

PYBIND11_MAKE_OPAQUE(std::vector<PerfectMatch>);
PYBIND11_MAKE_OPAQUE(std::vector<std::vector<PerfectMatch>>);

PYBIND11_MODULE(SDGpython, m) {
    py::enum_<NodeStatus>(m,"NodeStatus")
            .value("Active",NodeStatus::Active)
            .value("Deleted",NodeStatus::Deleted)
            .export_values();

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
        .def("__repr__",
             [](const LinkView &lv) {
                 return "<LinkView to " + std::to_string(lv.node().node_id()) + " at " + std::to_string(lv.distance()) + " bp>";
             })
         ;

    py::class_<NodeView>(m,"NodeView","A view for a Node in a Distance Graph")
        .def("node_id",&NodeView::node_id)
        .def("sequence",&NodeView::sequence)
        .def("size",&NodeView::size)
        .def("prev",&NodeView::prev,"A list with LinkViews to previous nodes",py::return_value_policy::move)
        .def("next",&NodeView::next,"A list with LinkViews to previous nodes",py::return_value_policy::move)
        .def("__repr__",
            [](const NodeView &nv) {
                return "< NodeView: " + std::to_string(nv.node_id()) + " on graph " + nv.graph().name + " bp>";
            })
        ;

    py::class_<DistanceGraph>(m, "DistanceGraph", "A Distance Graph")
        .def("get_all_nodeviews",&DistanceGraph::get_all_nodeviews,"include_disconnected"_a=true,"Returns a vector with NodeViews for all active nodes",py::return_value_policy::move)
        .def_readwrite("name",&DistanceGraph::name)
        ;

    py::class_<SequenceDistanceGraph,DistanceGraph>(m, "SequenceDistanceGraph", "A Sequence Distance Graph")
        ;

    py::class_<PairedReadsDatastore>(m, "PairedReadsDatastore", "A Paired Reads Datastore")
            .def("size",&PairedReadsDatastore::size)
            .def("get_read_sequence",&PairedReadsDatastore::get_read_sequence)
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

    py::class_<WorkSpace>(m, "WorkSpace", "A full SDG WorkSpace")
        .def(py::init<>(),"")
        .def(py::init<const std::string &>(),"Load from disk","filename"_a)
        .def_readonly("sdg",&WorkSpace::sdg)
        .def("add_long_reads_datastore",&WorkSpace::add_long_reads_datastore,"datastore"_a,"name"_a="",py::return_value_policy::reference)
        .def("add_paired_reads_datastore",&WorkSpace::add_paired_reads_datastore,"datastore"_a,"name"_a="",py::return_value_policy::reference)
        .def("add_linked_reads_datastore",&WorkSpace::add_linked_reads_datastore,"datastore"_a,"name"_a="",py::return_value_policy::reference)
        .def("get_long_reads_datastore",&WorkSpace::get_long_reads_datastore,"name"_a,py::return_value_policy::reference)
        .def("get_paired_reads_datastore",&WorkSpace::get_paired_reads_datastore,"name"_a,py::return_value_policy::reference)
        .def("get_linked_reads_datastore",&WorkSpace::get_linked_reads_datastore,"name"_a,py::return_value_policy::reference)
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
        .def_readonly("read_perfect_matches",&LongReadsRecruiter::read_perfect_matches)
        .def_readonly("node_reads",&LongReadsRecruiter::node_reads)
        .def("dump",&LongReadsRecruiter::dump,"filename"_a)
        .def("load",&LongReadsRecruiter::load,"filename"_a)
        .def("map",&LongReadsRecruiter::perfect_mappings,"hit_size"_a,"first_read"_a=1,"last_read"_a=0)
        ;


}

