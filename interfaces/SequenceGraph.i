%module bsg
%ignore Backbone;
%ignore PairedReadLinker;
%{
/* Includes the header in the wrapper code */
#include <cstdint>
#include <sstream>
#include "sglib/types/KmerTypes.hpp"
#include "sglib/types/GenericTypes.hpp"
#include "sglib/types/MappingTypes.hpp"
#include "sglib/graph/SequenceGraphPath.hpp"
#include "sglib/graph/SequenceSubGraph.hpp"
#include "sglib/graph/SequenceGraph.hpp"
#include "sglib/readers/SequenceGraphReader.h"
#include "sglib/processors/KmerCompressionIndex.hpp"
#include "sglib/workspace/WorkSpace.hpp"


#include "sglib/indexers/uniqueKmerIndex.hpp"
#include "sglib/mappers/LinkedReadMapper.hpp"
#include "sglib/mappers/LongReadMapper.hpp"
#include "sglib/mappers/PairedReadMapper.hpp"
#include "sglib/datastores/PairedReadsDatastore.hpp"
#include "sglib/datastores/PathsDatastore.hpp"

#include "sglib/mappers/threader/NodeMapper.h"
#include "sglib/mappers/threader/MappingThreader.h"

#include "sglib/processors/Untangler.hpp"
#include "sglib/processors/PathExplorer.h"
#include "sglib/graph/LinkageDiGraph.hpp"
#include "sglib/processors/LinkageUntangler.hpp"

%}


%include "python_docs.i"
%include "stdint.i"
%include "inttypes.i"
%include "std_string.i"
%include "std_vector.i"
%include "std_unordered_set.i"

%include "sglib/types/KmerTypes.hpp"
%include "sglib/types/GenericTypes.hpp"
%include "sglib/types/MappingTypes.hpp"
%include "sglib/graph/SequenceGraphPath.hpp"
%include "sglib/graph/SequenceSubGraph.hpp"
%include "sglib/graph/SequenceGraph.hpp"
%include "sglib/indexers/uniqueKmerIndex.hpp"
%include "sglib/processors/KmerCompressionIndex.hpp"
%include "sglib/workspace/WorkSpace.hpp"

%include "sglib/processors/Untangler.hpp"
%include "sglib/processors/PathExplorer.h"

%include "sglib/datastores/LinkedReadsDatastore.hpp"
%include "sglib/datastores/LongReadsDatastore.hpp"
%include "sglib/datastores/PairedReadsDatastore.hpp"
%include "sglib/datastores/PathsDatastore.hpp"
%include "sglib/mappers/LongReadMapper.hpp"
%include "sglib/mappers/LinkedReadMapper.hpp"
%include "sglib/mappers/PairedReadMapper.hpp"
%include "sglib/mappers/threader/NodeMapper.h"
%include "sglib/mappers/threader/MappingThreader.h"
%include "sglib/processors/LinkageUntangler.hpp"
%include "sglib/graph/LinkageDiGraph.hpp"


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
   %template(vectorUINT64) vector<uint64_t>;
   %template(vectorVectorUINT16) vector< vector <uint16_t > >;
   %template(vectorBool) vector<bool>;
   %template(vectorMinPosIDX) vector<MinPosIDX>;
   %template(vectorReadPosSize) vector< ReadPosSize >;
   %template(vectorReadMapping) vector<ReadMapping>;
   %template(vectorvectorReadMapping) vector<std::vector<ReadMapping>>;
   %template(vectorLongReadMapping) vector<LongReadMapping>;



   %ignore vector<LongReadMapper>::vector(size_type);
   %ignore vector<LongReadMapper>::resize;
   %template(vectorLongReadMapper) vector<LongReadMapper>;

   %ignore vector<LinkedReadMapper>::vector(size_type);
   %ignore vector<LinkedReadMapper>::resize;
   %template(vectorLinkedReadMapper) vector<LinkedReadMapper>;

   %ignore vector<PairedReadMapper>::vector(size_type);
   %ignore vector<PairedReadMapper>::resize;
   %template(vectorPairedReadMapper) vector<PairedReadMapper>;

   %ignore vector<PairedReadsDatastore>::vector(size_type);
   %ignore vector<PairedReadsDatastore>::resize;
   %template(vectorPairedReadsDatastore) vector<PairedReadsDatastore>;

   %ignore vector<PathsDatastore>::vector(size_type);
   %ignore vector<PathsDatastore>::resize;
   %template(vectorPathsDatastore) vector<PathsDatastore>;

   %template(SGNodePair) pair<sgNodeID_t, sgNodeID_t>;
   %template(vectorSGNodePair) vector<pair<sgNodeID_t, sgNodeID_t>>;
   %template(vectorvectorLink) vector<std::vector<Link>>;

   %ignore vector<SequenceGraphPath>::vector(size_type);
   %ignore vector<SequenceGraphPath>::resize;
   %template(vectorPath) vector<SequenceGraphPath>;
};

%feature("python:tp_hash") Link "std::hash<Link>";

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