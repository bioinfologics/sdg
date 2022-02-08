//
// Created by Luis Yanes (EI) on 22/03/2018.
//

#pragma once

#include <iostream>
#include <string>
#include <tuple>
#include <vector>
#include <limits>
#include <cstdint>
#include <xxhash/xxhash.h>
#include <set>
#include <unordered_set>
#include <algorithm>

using sgNodeID_t = int64_t; //first node is 1; negatives are RC
using seqID_t = int64_t; //first sequence is 0;

enum class NodeStatus:uint8_t {Active, Deleted};

enum class SupportType:uint8_t {Undefined,Operation,SequenceDistanceGraph,DistanceGraph,PairedRead,LinkedRead,LinkedTag,LongRead,ReadPath,KmerCoverage,StriderLink};

class Support{
public:
    Support(SupportType _type=SupportType::Undefined,uint16_t _index=0,int64_t _id=0):type(_type),index(_index),id(_id){};
    bool operator==(const Support &o) const {return std::tie(type,index,id)==std::tie(o.type,o.index,o.id);}
    SupportType type;
    uint16_t index;
    int64_t id;
};

/**
 * The Node contains the sequence of a node and its status {Active, Deleted}
 */
class Node{
public:
    Node(const std::string &_seq, NodeStatus _status) : sequence(_seq), status(_status){};
    Node(const std::string &_seq) : sequence(_seq),status(NodeStatus::Active){};
    Node() = default;
    bool operator==(const Node &o) const {
        return std::tie(status,sequence) == std::tie(o.status,o.sequence);
    }

    bool is_canonical();
    void make_rc();

    friend std::ostream &operator<<(std::ostream &os, const Node &node) {
        os << "Node ";
        if (node.status == NodeStatus::Deleted) os << "(deleted) ";
        os << node.sequence.size() <<" bp";
        return os;
    }

    std::string sequence = "";
    NodeStatus status = NodeStatus::Active;
    Support support;
};

/**
 * The Link represents a connection between one end of a sequence and the end of another sequence
 *
 * +AAAAAAAAA- connected to +BBBBBBBBBB-
 *
 * would be represented as the -A,B link. Links also contain a distance parameter in case of overlap the distance is
 * negative and in case of "scaffolding" link the distance would be positive.
 */
class Link{
public:
    Link(){};
    Link( sgNodeID_t _src, sgNodeID_t _dst, int32_t _dist, Support _support = {}) : source(_src), dest(_dst), dist(_dist), support(_support) {};

    bool operator==( const  Link) const;
    bool operator<(const Link)const;

    friend std::ostream &operator<<(std::ostream &os, const Link &link) {
        os << link.source << " -> " << link.dest;
        return os;
    }
    sgNodeID_t source = 0;
    sgNodeID_t dest = 0;
    int32_t dist = 0;
    Support support;
};


struct int128_hash {
    size_t operator()( const __int128 &x) const
    {
        const void *buffer = (unsigned char *) &x;
        uint64_t tmp_hash = XXH64(buffer, 16, 0);
        return tmp_hash;
    }
};

namespace sdglib {
    inline std::string str_rc(const std::string &sequence) {
        std::string rseq;
        rseq.resize(sequence.size());
        for (size_t i = 0, j = sequence.size() - 1; i < sequence.size(); ++i, --j) {
            switch (sequence[j]) {
                case 'A':
                    rseq[i] = 'T';
                    break;
                case 'C':
                    rseq[i] = 'G';
                    break;
                case 'G':
                    rseq[i] = 'C';
                    break;
                case 'T':
                    rseq[i] = 'A';
                    break;
                case 'a':
                    rseq[i] = 't';
                    break;
                case 'c':
                    rseq[i] = 'g';
                    break;
                case 'g':
                    rseq[i] = 'c';
                    break;
                case 't':
                    rseq[i] = 'a';
                    break;
            }
        }
        return rseq;
    }

    inline void reverse_path(std::vector<sgNodeID_t> &out,const std::vector<sgNodeID_t> &in){
        out.clear();
        out.reserve(in.size());
        for (auto it=in.crbegin();it<in.crend();++it) {
            out.emplace_back(-*it);
        }
    }

    inline std::vector<sgNodeID_t> reverse_path(const std::vector<sgNodeID_t> &in){
        std::vector<sgNodeID_t> out;
        sdglib::reverse_path(out,in);
        return out;
    }

    inline std::string kmer_to_sequence(__uint128_t kmer, uint8_t k) {
        std::string seq;
        seq.reserve(k);
        for (int shift=(k-1)*2;shift>=0;shift-=2) {
            //std::cout<<"kmer: "
            switch ((kmer>>shift)%4){
                case 0:
                    seq.push_back('A');
                    break;
                case 1:
                    seq.push_back('C');
                    break;
                case 2:
                    seq.push_back('G');
                    break;
                case 3:
                    seq.push_back('T');
                    break;
            }
        }
        return seq;
    }
}

namespace sdglib {
    template <class T>
            std::vector<T> vec_slice(const std::vector<T> &v, int64_t start=0, int64_t end=INT64_MAX){
                if (start<0) start+=v.size()-1;
                if (start<0) start=0;
                if (end<0) end+=v.size()-1;
                if (end<0) end=0;
                if (end>v.size()) end=v.size();
                if (start>=end) return {};
                return std::vector<T>(v.begin()+start,v.begin()+end);
            };
    template <class T>
            size_t intersection_size(const std::set<T>& v1, const std::set<T>& v2)
            {
                size_t s=0;
                for (auto p1=v1.cbegin(),p2=v2.cbegin();p1!=v1.cend() and p2!=v2.cend();){
                    if (*p1==*p2) {
                        ++s;
                        ++p1;
                        ++p2;
                    }
                    else if (*p1<*p2) ++p1;
                    else ++p2;
                }
                return s;
            };
    template <class T>
            size_t intersection_size(const std::unordered_set<T>& v1, const std::unordered_set<T>& v2)
            {
                size_t s=0;
                for (const auto & x:v1) if (v2.count(x)) ++s;
                return s;
            };
    template <class T>
            T rough_median(std::vector<T> v)
            {
                if (v.empty()) throw std::runtime_error("tried to compute rough_median of empty vector");
                std::sort(v.begin(),v.end());
                return v[v.size()/2];
            }
}