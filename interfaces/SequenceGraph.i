%module bsg
%ignore Backbone;
%ignore PairedReadLinker;
%{
/* Includes the header in the wrapper code */
#include <cstdint>
#include <sstream>
#include "sglib/types/KmerTypes.hpp"
#include "sglib/types/GenericTypes.hpp"
#include "sglib/graph/SequenceGraphPath.hpp"
#include "sglib/graph/SequenceSubGraph.hpp"
#include "sglib/graph/SequenceGraph.hpp"
#include "sglib/readers/SequenceGraphReader.h"
#include "sglib/KmerCompressionIndex.hpp"
#include "sglib/workspace/WorkSpace.hpp"


#include "sglib/indexers/uniqueKmerIndex.hpp"
#include "sglib/mappers/LinkedReadMapper.hpp"
#include "sglib/mappers/LongReadMapper.hpp"
#include "sglib/mappers/PairedReadMapper.hpp"
#include "sglib/mappers/threader/NodeMapper.h"
#include "sglib/mappers/threader/MappingThreader.h"

#include "sglib/processors/Untangler.hpp"
%}

%include "python_docs.i"
%include "stdint.i"
%include "inttypes.i"
%include "std_string.i"
%include "std_vector.i"
%include "std_unordered_set.i"
%include "sglib/types/KmerTypes.hpp"
%include "sglib/types/GenericTypes.hpp"
%include "sglib/graph/SequenceGraphPath.hpp"
%include "sglib/graph/SequenceSubGraph.hpp"
%include "sglib/graph/SequenceGraph.hpp"
%include "sglib/indexers/uniqueKmerIndex.hpp"
%include "sglib/KmerCompressionIndex.hpp"
%include "sglib/workspace/WorkSpace.hpp"

%include "sglib/processors/Untangler.hpp"

%include "sglib/datastores/LinkedReadsDatastore.hpp"
%include "sglib/datastores/LongReadsDatastore.hpp"
%include "sglib/datastores/PairedReadsDatastore.hpp"
%include "sglib/mappers/LongReadMapper.hpp"
%include "sglib/mappers/LinkedReadMapper.hpp"
%include "sglib/mappers/PairedReadMapper.hpp"
%include "sglib/mappers/threader/NodeMapper.h"
%include "sglib/mappers/threader/MappingThreader.h"

namespace std {
   %template(vectorInt) vector<int>;
   %template(vectorDouble) vector<double>;
   %template(vectorFloat) vector<float>;
   %template(vectorString) vector<string>;
   %template(vectorLink) vector<Link>;
   %template(vectorSGNode) vector<sgNodeID_t>;
   %template(vectorNode) vector<Node>;
   %template(vectorNodeVisitor) vector<nodeVisitor>;
   %template(vectorUINT16) vector<uint16_t >;
   %template(vectorVectorUINT16) vector< vector <uint16_t > >;
   %template(vectorBool) vector<bool>;
   %template(vectorMinPosIDX) vector<MinPosIDX>;
   %template(vectorReadPosSize) std::vector< ReadPosSize >;
   %template(vectorReadMapping) std::vector<ReadMapping>;
   %template(vectorLongReadMapping) std::vector<LongReadMapping>;
   %template(SGNodePair) std::pair<sgNodeID_t, sgNodeID_t>;
   %template(vectorSGNodePair) std::vector<std::pair<sgNodeID_t, sgNodeID_t>>;
   %template(vectorvectorLink) std::vector<std::vector<Link>>;
   %ignore vector<SequenceGraphPath>::vector(size_type);
   %ignore vector<SequenceGraphPath>::resize;
   %template(vectorPath) vector<SequenceGraphPath>;
};


%define __STR__()
std::string __str__() {
  std::ostringstream out;
  out << *$self;
  return out.str().c_str();
}
%enddef

%extend Link{
    __STR__();
}

%extend Node{
    __STR__();
}

//TODO: Make sure this changes on new releases
%pythoncode %{
__version__ = "0.1"
%}