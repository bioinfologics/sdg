%module bsg
%{
/* Includes the header in the wrapper code */
#include <cstdint>
#include "sglib/types/KmerTypes.hpp"
#include "sglib/types/GenericTypes.hpp"
#include "sglib/graph/SequenceSubGraph.hpp"
#include "sglib/graph/SequenceGraphPath.hpp"
#include "sglib/graph/SequenceGraph.hpp"
#include "sglib/readers/SequenceGraphReader.h"
#include "sglib/KmerCompressionIndex.hpp"
#include "sglib/WorkSpace.hpp"
%}
%include "stdint.i"
%include "inttypes.i"
%include "std_string.i"
%include "std_vector.i"
%include "std_unordered_set.i"
%include "sglib/types/KmerTypes.hpp"
%include "sglib/types/GenericTypes.hpp"
%include "sglib/graph/SequenceSubGraph.hpp"
%include "sglib/graph/SequenceGraphPath.hpp"
/* Parse the header file to generate wrappers */
%include "sglib/graph/SequenceGraph.hpp"
%include "sglib/KmerCompressionIndex.hpp"
%include "sglib/WorkSpace.hpp"

namespace std {
   %template(vectorInt) vector<int>;
   %template(vectorDouble) vector<double>;
   %template(vectorFloat) vector<float>;
   %template(vectorString) vector<string>;
   %template(vectorSGNode) vector<sgNodeID_t>;
//   %template(vectorPath) vector<SequenceGraphPath>;
   %template(vectorNode) vector<Node>;
   %template(vectorNodeVisitor) vector<nodeVisitor>;
   %template(vectorUINT16) vector<uint16_t >;
   %template(vectorVectorUINT16) vector< vector <uint16_t > >;
   %template(vectorBool) vector<bool>;
};

