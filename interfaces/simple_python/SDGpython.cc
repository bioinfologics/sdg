//
// Created by Bernardo Clavijo (EI) on 13/02/2020.
//

#include <deps/pybind11/pybind11.h>
#include <deps/pybind11/stl.h>
#include <sdglib/graph/SequenceDistanceGraph.hpp>
#include <sdglib/workspace/WorkSpace.hpp>
#include <sdglib/views/NodeView.hpp>

namespace py = pybind11;
using namespace py::literals;



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
        .def("next",&NodeView::next,"A list with LinkViews to previous nodes",py::return_value_policy::move);

    py::class_<SequenceDistanceGraph>(m, "SequenceDistanceGraph", "A Sequence Distance Graph")
        .def("get_all_nodeviews",&SequenceDistanceGraph::get_all_nodeviews,"Returns a vector with NodeViews for all active nodes",py::return_value_policy::move)
        //.def("name",&SequenceDistanceGraph::name)
        ;

    py::class_<WorkSpace>(m, "WorkSpace", "A full SDG WorkSpace")
        .def(py::init<>(),"")
        .def(py::init<const std::string &>(),"Load from disk","filename"_a)
        .def_readonly("sdg",&WorkSpace::sdg)
        ;
}