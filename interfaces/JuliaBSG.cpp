//
// Created by Ben Ward (EI) on 02/10/2018.
//

#include "jlcxx/jlcxx.hpp"
#include <sglib/graph/SequenceGraph.hpp>

JLCXX_MODULE define_julia_module(jlcxx::Module& mod)
{
    mod.add_type<SequenceGraph>("SequenceGraph").method("load_from_gfa", &SequenceGraph::load_from_gfa);

}