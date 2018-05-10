%module bsg
%{
/* Includes the header in the wrapper code */
#include "sglib/types/KmerTypes.hpp"
#include "sglib/types/GenericTypes.hpp"
#include "sglib/graph/SequenceSubGraph.hpp"
#include "sglib/graph/SequenceGraphPath.hpp"
#include "sglib/graph/SequenceGraph.hpp"
%}
%include "stdint.i"
%include "std_string.i"
%include "std_vector.i"
%include "std_unordered_set.i"
%include "sglib/types/KmerTypes.hpp"
%include "sglib/types/GenericTypes.hpp"
%include "sglib/graph/SequenceSubGraph.hpp"
%include "sglib/graph/SequenceGraphPath.hpp"
/* Parse the header file to generate wrappers */
%include "sglib/graph/SequenceGraph.hpp"

namespace std {
   %template(vectorInt) vector<int>;
   %template(vectorDouble) vector<double>;
   %template(vectorFloat) vector<float>;
   %template(vectorString) vector<string>;
   %template(vectorSGNode) vector<sgNodeID_t>;
//   %template(vectorPath) vector<SequenceGraphPath>;
   %template(vectorNode) vector<Node>;
   %template(vectorNodeVisitor) vector<nodeVisitor>;
};

