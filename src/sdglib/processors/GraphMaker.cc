//
// Created by Bernardo Clavijo (EI) on 05/07/2018.
//

#include <sdglib/bloom/BloomFilter.hpp>
#include <sdglib/batch_counter/BatchKmersCounter.hpp>
#include "GraphMaker.hpp"

struct sbcont {
    unsigned int fA:1;
    unsigned int fC:1;
    unsigned int fG:1;
    unsigned int fT:1;
    unsigned int bA:1;
    unsigned int bC:1;
    unsigned int bG:1;
    unsigned int bT:1;
};

struct fb {
    unsigned int fw:4;
    unsigned int bw:4;
};
union kcontext {
    struct fb fb;
    struct sbcont sbcont;
};

std::string kmer_to_sequence(uint64_t kmer, uint8_t k) {
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

std::array<uint64_t,4> kmer_fw_neighbours(uint64_t kmer, uint8_t k) {
    std::array<uint64_t,4> n;
    for (uint64_t i=0;i<4;++i){
        n[i]=((kmer<<2)+i)%(((uint64_t) 1)<<(k*2));
    }
//    std::cout<<"Produced FW neighbours for "<<kmer_to_sequence(kmer,k)<<std::endl;
//    std::cout<<" for A: "<<kmer_to_sequence(n[0],k)<<std::endl;
//    std::cout<<" for C: "<<kmer_to_sequence(n[1],k)<<std::endl;
//    std::cout<<" for G: "<<kmer_to_sequence(n[2],k)<<std::endl;
//    std::cout<<" for T: "<<kmer_to_sequence(n[3],k)<<std::endl;
    return n;
}

std::array<uint64_t,4> kmer_bw_neighbours(uint64_t kmer, uint8_t k) {
    std::array<uint64_t,4> n;
    for (uint64_t i=0;i<4;++i){
        n[i]=(kmer>>2)+(((uint64_t) i)<<((k-1)*2));
    }
//    std::cout<<"Produced BW neighbours for "<<kmer_to_sequence(kmer,k)<<std::endl;
//    std::cout<<" for A: "<<kmer_to_sequence(n[0],k)<<std::endl;
//    std::cout<<" for C: "<<kmer_to_sequence(n[1],k)<<std::endl;
//    std::cout<<" for G: "<<kmer_to_sequence(n[2],k)<<std::endl;
//    std::cout<<" for T: "<<kmer_to_sequence(n[3],k)<<std::endl;
    return n;
}

uint64_t kmer_reverse(uint64_t kmer, uint8_t k) {
    uint64_t old_kmer=kmer;
    uint64_t rkmer=0;
    for (int pos=0;pos<k;++pos) {
        //std::cout<<"kmer: "
        //3->0 T->A     2->1 G C     1->2
        uint64_t nuc=3-(old_kmer%4);
        old_kmer>>=2;
        //std::cout<<"rc adding "<<nuc<<std::endl;
        rkmer<<=2;
        rkmer+=nuc;
    }
    return rkmer;
}

uint64_t kmer_cannonical(uint64_t kmer, uint8_t k) {
    uint64_t rkmer=kmer_reverse(kmer,k);
    return (rkmer<kmer ? rkmer:kmer);
}

uint64_t kmer_back_ovl(uint64_t kmer, uint8_t k){
    return kmer>>2;
}

uint64_t kmer_fw_ovl(uint64_t kmer, uint8_t k){
    return  kmer%(((uint64_t) 1)<<((k-1)*2));
}


std::string kmer_to_sequence128(__uint128_t kmer, uint8_t k) {
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

std::array<__uint128_t,4> kmer_fw_neighbours128(__uint128_t kmer, uint8_t k) {
    std::array<__uint128_t,4> n;
    for (__uint128_t i=0;i<4;++i){
        n[i]=((kmer<<2)+i)%(((__uint128_t) 1)<<(k*2));
    }
    return n;
}

std::array<__uint128_t,4> kmer_bw_neighbours128(__uint128_t kmer, uint8_t k) {
    std::array<__uint128_t,4> n;
    for (__uint128_t i=0;i<4;++i){
        n[i]=(kmer>>2)+(((__uint128_t) i)<<((k-1)*2));
    }
    return n;
}

__uint128_t kmer_reverse128(__uint128_t kmer, uint8_t k) {
    __uint128_t old_kmer=kmer;
    __uint128_t rkmer=0;
    for (int pos=0;pos<k;++pos) {
        __uint128_t nuc=3-(old_kmer%4);
        old_kmer>>=2;
        rkmer<<=2;
        rkmer+=nuc;
    }
    return rkmer;
}

__uint128_t kmer_cannonical128(__uint128_t kmer, uint8_t k) {
    __uint128_t rkmer=kmer_reverse128(kmer,k);
    return (rkmer<kmer ? rkmer:kmer);
}

class klidxs128 {
public:
    klidxs128(__uint128_t _kmer, uint64_t _idx): kmer(_kmer), idx(_idx){};
    __uint128_t kmer;
    uint64_t idx;
};

inline void get_fw_idxs128(std::vector<klidxs128> &out, const __uint128_t kmer, uint8_t k, const std::vector<__uint128_t> & kmerlist){
    out.clear();
    __uint128_t cnext;
    //std::vector<__uint128_t>::const_iterator cnitr;
    for (const __uint128_t &fnk:kmer_fw_neighbours128(kmer, k)) {
        cnext=kmer_cannonical128(fnk,k);
        auto cnitr=std::lower_bound(kmerlist.begin(),kmerlist.end(),cnext);
        if (*cnitr==cnext)
            out.emplace_back(fnk,cnitr-kmerlist.begin());
    }
}

class klidxs64 {
public:
    klidxs64(uint64_t _kmer, uint64_t _idx): kmer(_kmer), idx(_idx){};
    uint64_t kmer;
    uint64_t idx;
};

inline void get_fw_idxs64(std::vector<klidxs64> &out, const uint64_t kmer, uint8_t k, const std::vector<uint64_t> & kmerlist){
    out.clear();
    uint64_t cnext;
    //std::vector<__uint128_t>::const_iterator cnitr;
    for (const uint64_t &fnk:kmer_fw_neighbours(kmer, k)) {
        cnext=kmer_cannonical(fnk,k);
        auto cnitr=std::lower_bound(kmerlist.begin(),kmerlist.end(),cnext);
        if (*cnitr==cnext)
            out.emplace_back(fnk,cnitr-kmerlist.begin());
    }
}

template<class T>
inline bool in_sorted_vector(const std::vector<T> &V,T VAL) {
    auto itr=lower_bound(V.begin(),V.end(),VAL);
    return itr!=V.end() and *itr==VAL;
}

void GraphMaker::new_graph_from_kmerlist_128(const std::vector<__uint128_t> & kmerlist,uint8_t k) {
    //TODO: add a bloom filter to speed this up
    sg.nodes.clear();
    sg.links.clear();
    sg.oldnames.clear();
    sg.add_node(Node("",NodeStatus::Deleted));
    sdglib::OutputLog()<<"Constructing Graph from "<<kmerlist.size()<<" "<<std::to_string(k)<<"-mers"<<std::endl;
    std::string s;
    s.reserve(1000000); //avoid contig-sequence growth
    std::vector<bool> used_kmers(kmerlist.size());
    std::vector<kcontext> kcontextlist(kmerlist.size());
    std::vector<bool> is_end_fw(kmerlist.size());
    std::vector<bool> is_end_bw(kmerlist.size());
    sdglib::OutputLog()<<"Finding neighbours"<<std::endl;
#pragma omp parallel for
    for (uint64_t i=0;i<kmerlist.size();++i){
        auto fns=kmer_fw_neighbours128(kmerlist[i], k);
        auto &c=kcontextlist[i];
        c.fb.fw=0;
        c.fb.bw=0;
        if (in_sorted_vector(kmerlist,kmer_cannonical128(fns[0],k))) c.sbcont.fA=1;
        if (in_sorted_vector(kmerlist,kmer_cannonical128(fns[1],k))) c.sbcont.fC=1;
        if (in_sorted_vector(kmerlist,kmer_cannonical128(fns[2],k))) c.sbcont.fG=1;
        if (in_sorted_vector(kmerlist,kmer_cannonical128(fns[3],k))) c.sbcont.fT=1;
        auto bns=kmer_bw_neighbours128(kmerlist[i], k);
        if (in_sorted_vector(kmerlist,kmer_cannonical128(bns[0],k))) c.sbcont.bA=1;
        if (in_sorted_vector(kmerlist,kmer_cannonical128(bns[1],k))) c.sbcont.bC=1;
        if (in_sorted_vector(kmerlist,kmer_cannonical128(bns[2],k))) c.sbcont.bG=1;
        if (in_sorted_vector(kmerlist,kmer_cannonical128(bns[3],k))) c.sbcont.bT=1;
    }
    sdglib::OutputLog()<<"Marking ends"<<std::endl;
#pragma omp parallel for
    for (uint64_t i=0;i<kmerlist.size();++i) {
        auto &c=kcontextlist[i];
        if (c.fb.bw!=1 and c.fb.bw!=2 and c.fb.bw!=4 and c.fb.bw!=8 ) {
            is_end_bw[i]=true;
        } else {
            int nn;
            if (c.sbcont.bA) nn=0;
            else if (c.sbcont.bC) nn=1;
            else if (c.sbcont.bG) nn=2;
            else if (c.sbcont.bT) nn=3;
            auto next=kmer_bw_neighbours128(kmerlist[i],k)[nn];
            auto cnext=kmer_cannonical128(next,k);
            auto ni=std::lower_bound(kmerlist.begin(),kmerlist.end(),cnext)-kmerlist.begin();
            auto &nc=kcontextlist[ni];
            if (next==cnext and nc.fb.fw!=1 and nc.fb.fw!=2 and nc.fb.fw!=4 and nc.fb.fw!=8) is_end_bw[i]=true;
            if (next!=cnext and nc.fb.bw!=1 and nc.fb.bw!=2 and nc.fb.bw!=4 and nc.fb.bw!=8) is_end_bw[i]=true;
        }
        if (c.fb.fw!=1 and c.fb.fw!=2 and c.fb.fw!=4 and c.fb.fw!=8 ) {
            is_end_fw[i]=true;
        } else {
            int nn;
            if (c.sbcont.fA) nn=0;
            else if (c.sbcont.fC) nn=1;
            else if (c.sbcont.fG) nn=2;
            else if (c.sbcont.fT) nn=3;
            auto next=kmer_fw_neighbours128(kmerlist[i],k)[nn];
            auto cnext=kmer_cannonical128(next,k);
            auto ni=std::lower_bound(kmerlist.begin(),kmerlist.end(),cnext)-kmerlist.begin();
            auto &nc=kcontextlist[ni];
            if (next==cnext and nc.fb.bw!=1 and nc.fb.bw!=2 and nc.fb.bw!=4 and nc.fb.bw!=8) is_end_fw[i]=true;
            if (next!=cnext and nc.fb.fw!=1 and nc.fb.fw!=2 and nc.fb.fw!=4 and nc.fb.fw!=8) is_end_fw[i]=true;
        }
    }
    sdglib::OutputLog()<<"Creating unitigs"<<std::endl;

    for (uint64_t start_kmer_idx=0;start_kmer_idx<kmerlist.size();++start_kmer_idx) {
        if (used_kmers[start_kmer_idx]) continue; //any kmer can only belong to one unitig.ww

        //Check this k-mer is an end/junction
        auto start_kmer = kmerlist[start_kmer_idx];


        auto end_bw=is_end_bw[start_kmer_idx];
        auto end_fw=is_end_fw[start_kmer_idx];

        if (!end_bw and !end_fw) continue;

        if (end_bw and end_fw) {
            //kmer as unitig
            s=kmer_to_sequence128(start_kmer, k);
            used_kmers[start_kmer_idx] = true;
        } else {
            //unitig start on kmer
            //make sure the unitig starts on FW

            auto current_kmer = start_kmer;
            auto current_idx = start_kmer_idx;
            used_kmers[start_kmer_idx] = true;

            if (end_fw) {
                current_kmer = kmer_reverse128(start_kmer, k);
                std::swap(end_fw,end_bw);
            }

            //Add kmer as unitig
            s=kmer_to_sequence128(current_kmer, k);

            unsigned char nucleotides[4] = {'A', 'C', 'G', 'T'};



            //std::cout<<"Starting sequence construction at kmer "<<kmer_to_sequence128(current_kmer,k)<<std::endl;
            while (!end_fw) {
                //Add end nucleotide, update current_kmer;
                //get_fw_idxs128(fwn, current_kmer, k, kmerlist);
                auto c=kcontextlist[current_idx];
                int nn;
                if (kmer_cannonical128(current_kmer,k)==current_kmer) {
                    if (c.sbcont.fA) nn = 0;
                    else if (c.sbcont.fC) nn = 1;
                    else if (c.sbcont.fG) nn = 2;
                    else if (c.sbcont.fT) nn = 3;
                }
                else{
                    if (c.sbcont.bA) nn = 3;
                    else if (c.sbcont.bC) nn = 2;
                    else if (c.sbcont.bG) nn = 1;
                    else if (c.sbcont.bT) nn = 0;
                }
                current_kmer=kmer_fw_neighbours128(current_kmer,k)[nn];
                auto ccurr=kmer_cannonical128(current_kmer,k);
                current_idx=std::lower_bound(kmerlist.begin(),kmerlist.end(),ccurr)-kmerlist.begin();
                if (used_kmers[current_idx]) break; //this should never happen, since these have ends.
                used_kmers[current_idx]=true;
                s.push_back(nucleotides[current_kmer % 4]);
                if (ccurr==current_kmer) end_fw = is_end_fw[current_idx];
                else end_fw = is_end_bw[current_idx];
            }
        }
        sg.add_node(Node(s));
        if (!sg.nodes.back().is_canonical()) sg.nodes.back().make_rc();//inefficient, but once in a node
    }
    sdglib::OutputLog()<<sg.nodes.size()<<" unitigs"<<std::endl;
    //If there are any perfect circles, they won't have ends, so just pick any unused kmer and go fw till you find the same k-mer fw.
    //(this is going to be even tricker to parallelise)
    sdglib::OutputLog()<<"doing the circles now"<<std::endl;
    for (uint64_t start_kmer_idx=0;start_kmer_idx<kmerlist.size();++start_kmer_idx) {
        if (used_kmers[start_kmer_idx]) continue; //any kmer can only belong to one unitig.ww

        //Check this k-mer is an end/junction
        auto start_kmer = kmerlist[start_kmer_idx];
        std::vector<klidxs128> fwn;
        unsigned char nucleotides[4] = {'A', 'C', 'G', 'T'};
        auto current_kmer = start_kmer;
        used_kmers[start_kmer_idx] = true;
        s=kmer_to_sequence128(current_kmer, k);

        while (true) {
            //Add end nucleotide, update current_kmer;
            get_fw_idxs128(fwn, current_kmer, k, kmerlist);
            current_kmer = fwn[0].kmer;
            if (current_kmer==start_kmer) break;
            used_kmers[fwn[0].idx] = true;
            s.push_back(nucleotides[current_kmer % 4]);
        }
        sg.add_node(Node(s));
        if (!sg.nodes.back().is_canonical()) sg.nodes.back().make_rc();//inefficient, but once in a node
    }
    sdglib::OutputLog()<<sg.nodes.size()<<" unitigs"<<std::endl;

    //save the (k-1)mer in (rev on first k-1 / fw on last k-1) or out ( fw on first k-1 / bw on last k-1)
    std::vector<std::pair<__uint128_t,sgNodeID_t>> in, out;
    in.reserve(2*sg.nodes.size());
    out.reserve(2*sg.nodes.size());
    CStringKMerFactory128 skf_ovl(k-1);
    for (uint64_t nid=1;nid<sg.nodes.size();++nid){
        std::vector<std::pair<__uint128_t,bool>> first,last;
        skf_ovl.create_kmers_direction(first,sg.nodes[nid].sequence.substr(0,k-1).c_str());
        if (first[0].second) out.emplace_back(first[0].first,nid);
        else in.emplace_back(first[0].first,nid);
        skf_ovl.create_kmers_direction(last,sg.nodes[nid].sequence.substr(sg.nodes[nid].sequence.size()-k+1,k-1).c_str());
        if (last[0].second) in.emplace_back(last[0].first,-nid);
        else out.emplace_back(last[0].first,-nid);
    }
    std::sort(in.begin(),in.end());
    std::sort(out.begin(),out.end());
    //connect out->in for all combinations on each kmer
    uint64_t next_out_idx=0;
    for(auto &i:in){
        while(next_out_idx<out.size() and out[next_out_idx].first<i.first) ++next_out_idx;
        for (auto oidx=next_out_idx;oidx<out.size() and out[oidx].first==i.first;++oidx) {
            sg.add_link(i.second,out[oidx].second,-k+1); //no support, although we could add the DBG operation as such
        }
    }
    sdglib::OutputLog()<<"Graph construction finished"<<std::endl;
}

void GraphMaker::new_graph_from_kmerlist_64(const std::vector<__uint64_t> & kmerlist,uint8_t k) {
    //TODO: add a bloom filter to speed this up
    sg.nodes.clear();
    sg.links.clear();
    sg.oldnames.clear();
    sg.add_node(Node("",NodeStatus::Deleted));
    sdglib::OutputLog()<<"Constructing Graph from "<<kmerlist.size()<<" "<<std::to_string(k)<<"-mers"<<std::endl;
    std::string s;
    s.reserve(1000000); //avoid contig-sequence growth
    std::vector<bool> used_kmers(kmerlist.size());
    std::vector<kcontext> kcontextlist(kmerlist.size());
    std::vector<bool> is_end_fw(kmerlist.size());
    std::vector<bool> is_end_bw(kmerlist.size());
    sdglib::OutputLog()<<"Finding neighbours"<<std::endl;
#pragma omp parallel for
    for (uint64_t i=0;i<kmerlist.size();++i){
        auto fns=kmer_fw_neighbours(kmerlist[i], k);
        auto &c=kcontextlist[i];
        c.fb.fw=0;
        c.fb.bw=0;
        if (in_sorted_vector(kmerlist,kmer_cannonical(fns[0],k))) c.sbcont.fA=1;
        if (in_sorted_vector(kmerlist,kmer_cannonical(fns[1],k))) c.sbcont.fC=1;
        if (in_sorted_vector(kmerlist,kmer_cannonical(fns[2],k))) c.sbcont.fG=1;
        if (in_sorted_vector(kmerlist,kmer_cannonical(fns[3],k))) c.sbcont.fT=1;
        auto bns=kmer_bw_neighbours(kmerlist[i], k);
        if (in_sorted_vector(kmerlist,kmer_cannonical(bns[0],k))) c.sbcont.bA=1;
        if (in_sorted_vector(kmerlist,kmer_cannonical(bns[1],k))) c.sbcont.bC=1;
        if (in_sorted_vector(kmerlist,kmer_cannonical(bns[2],k))) c.sbcont.bG=1;
        if (in_sorted_vector(kmerlist,kmer_cannonical(bns[3],k))) c.sbcont.bT=1;
    }
    sdglib::OutputLog()<<"Marking ends"<<std::endl;
#pragma omp parallel for
    for (uint64_t i=0;i<kmerlist.size();++i) {
        auto &c=kcontextlist[i];
        if (c.fb.bw!=1 and c.fb.bw!=2 and c.fb.bw!=4 and c.fb.bw!=8 ) {
            is_end_bw[i]=true;
        } else {
            int nn;
            if (c.sbcont.bA) nn=0;
            else if (c.sbcont.bC) nn=1;
            else if (c.sbcont.bG) nn=2;
            else if (c.sbcont.bT) nn=3;
            auto next=kmer_bw_neighbours(kmerlist[i],k)[nn];
            auto cnext=kmer_cannonical(next,k);
            auto ni=std::lower_bound(kmerlist.begin(),kmerlist.end(),cnext)-kmerlist.begin();
            auto &nc=kcontextlist[ni];
            if (next==cnext and nc.fb.fw!=1 and nc.fb.fw!=2 and nc.fb.fw!=4 and nc.fb.fw!=8) is_end_bw[i]=true;
            if (next!=cnext and nc.fb.bw!=1 and nc.fb.bw!=2 and nc.fb.bw!=4 and nc.fb.bw!=8) is_end_bw[i]=true;
        }
        if (c.fb.fw!=1 and c.fb.fw!=2 and c.fb.fw!=4 and c.fb.fw!=8 ) {
            is_end_fw[i]=true;
        } else {
            int nn;
            if (c.sbcont.fA) nn=0;
            else if (c.sbcont.fC) nn=1;
            else if (c.sbcont.fG) nn=2;
            else if (c.sbcont.fT) nn=3;
            auto next=kmer_fw_neighbours(kmerlist[i],k)[nn];
            auto cnext=kmer_cannonical(next,k);
            auto ni=std::lower_bound(kmerlist.begin(),kmerlist.end(),cnext)-kmerlist.begin();
            auto &nc=kcontextlist[ni];
            if (next==cnext and nc.fb.bw!=1 and nc.fb.bw!=2 and nc.fb.bw!=4 and nc.fb.bw!=8) is_end_fw[i]=true;
            if (next!=cnext and nc.fb.fw!=1 and nc.fb.fw!=2 and nc.fb.fw!=4 and nc.fb.fw!=8) is_end_fw[i]=true;
        }
    }
    uint64_t eb=0;
    for (const auto &e:is_end_bw) if(e) ++eb;
    uint64_t ef=0;
    for (const auto &e:is_end_fw) if(e) ++ef;
    sdglib::OutputLog()<<is_end_bw.size()<<" k-mers, with "<<eb<<" ends backward and "<<ef <<" ends forward"<<std::endl;
    sdglib::OutputLog()<<"Creating unitigs"<<std::endl;

    for (uint64_t start_kmer_idx=0;start_kmer_idx<kmerlist.size();++start_kmer_idx) {
        if (used_kmers[start_kmer_idx]) continue; //any kmer can only belong to one unitig.ww

        //Check this k-mer is an end/junction
        auto start_kmer = kmerlist[start_kmer_idx];


        auto end_bw=is_end_bw[start_kmer_idx];
        auto end_fw=is_end_fw[start_kmer_idx];

        if (!end_bw and !end_fw) continue;

        if (end_bw and end_fw) {
            //kmer as unitig
            s=kmer_to_sequence(start_kmer, k);
            used_kmers[start_kmer_idx] = true;
        } else {
            //unitig start on kmer
            //make sure the unitig starts on FW

            auto current_kmer = start_kmer;
            auto current_idx = start_kmer_idx;
            used_kmers[start_kmer_idx] = true;

            if (end_fw) {
                current_kmer = kmer_reverse(start_kmer, k);
                std::swap(end_fw,end_bw);
            }

            //Add kmer as unitig
            s=kmer_to_sequence(current_kmer, k);

            unsigned char nucleotides[4] = {'A', 'C', 'G', 'T'};



            //std::cout<<"Starting sequence construction at kmer "<<kmer_to_sequence128(current_kmer,k)<<std::endl;
            while (!end_fw) {
                //Add end nucleotide, update current_kmer;
                //get_fw_idxs128(fwn, current_kmer, k, kmerlist);
                auto c=kcontextlist[current_idx];
                int nn;
                if (kmer_cannonical(current_kmer,k)==current_kmer) {
                    if (c.sbcont.fA) nn = 0;
                    else if (c.sbcont.fC) nn = 1;
                    else if (c.sbcont.fG) nn = 2;
                    else if (c.sbcont.fT) nn = 3;
                }
                else{
                    if (c.sbcont.bA) nn = 3;
                    else if (c.sbcont.bC) nn = 2;
                    else if (c.sbcont.bG) nn = 1;
                    else if (c.sbcont.bT) nn = 0;
                }
                current_kmer=kmer_fw_neighbours(current_kmer,k)[nn];
                auto ccurr=kmer_cannonical(current_kmer,k);
                current_idx=std::lower_bound(kmerlist.begin(),kmerlist.end(),ccurr)-kmerlist.begin();
                if (used_kmers[current_idx]) break; //this should never happen, since these have ends.
                used_kmers[current_idx]=true;
                s.push_back(nucleotides[current_kmer % 4]);
                if (ccurr==current_kmer) end_fw = is_end_fw[current_idx];
                else end_fw = is_end_bw[current_idx];
            }
        }
        sg.add_node(Node(s));
        if (!sg.nodes.back().is_canonical()) sg.nodes.back().make_rc();//inefficient, but once in a node
    }
    sdglib::OutputLog()<<sg.nodes.size()<<" unitigs"<<std::endl;
    uint64_t uk=0;
    for (auto const & u:used_kmers) if(u) ++uk;
    sdglib::OutputLog()<<uk<<"/"<<used_kmers.size()<<" kmers are used"<<std::endl;
    //If there are any perfect circles, they won't have ends, so just pick any unused kmer and go fw till you find the same k-mer fw.
    //(this is going to be even tricker to parallelise)
    sdglib::OutputLog()<<"doing the circles now"<<std::endl;
    uint64_t cid=0;
    for (uint64_t start_kmer_idx=0;start_kmer_idx<kmerlist.size();++start_kmer_idx) {
        if (used_kmers[start_kmer_idx]) continue; //any kmer can only belong to one unitig.ww
        std::cout<<"Circle #"<<cid++<<"started on kmer "<<start_kmer_idx<<": "<<kmer_to_sequence(kmerlist[start_kmer_idx],k)<<std::endl;
        //Check this k-mer is an end/junction
        auto start_kmer = kmerlist[start_kmer_idx];
        std::vector<klidxs64> fwn;
        unsigned char nucleotides[4] = {'A', 'C', 'G', 'T'};
        auto current_kmer = start_kmer;
        used_kmers[start_kmer_idx] = true;
        s=kmer_to_sequence(current_kmer, k);

        while (true) {
            //Add end nucleotide, update current_kmer;
            get_fw_idxs64(fwn, current_kmer, k, kmerlist);
            current_kmer = fwn[0].kmer;
            std::cout<<"Circle goes on to kmer "<<fwn[0].idx<<": "<<kmer_to_sequence(fwn[0].kmer,k)<<std::endl;
            if (current_kmer==start_kmer) break;
            if (used_kmers[fwn[0].idx]){
                std::cout<<"ERROR: using a kmer that is already used, but i'm supposed to be on a circle"<<std::endl;
                break;
            }
            used_kmers[fwn[0].idx] = true;
            s.push_back(nucleotides[current_kmer % 4]);
        }
        sg.add_node(Node(s));
        if (!sg.nodes.back().is_canonical()) sg.nodes.back().make_rc();//inefficient, but once in a node
    }
    sdglib::OutputLog()<<sg.nodes.size()<<" unitigs"<<std::endl;

    //save the (k-1)mer in (rev on first k-1 / fw on last k-1) or out ( fw on first k-1 / bw on last k-1)
    std::vector<std::pair<uint64_t,sgNodeID_t>> in, out;
    in.reserve(2*sg.nodes.size());
    out.reserve(2*sg.nodes.size());
    CStringKMerFactory skf_ovl(k-1);
    for (uint64_t nid=1;nid<sg.nodes.size();++nid){
        std::vector<std::pair<uint64_t,bool>> first,last;
        skf_ovl.create_kmers_direction(first,sg.nodes[nid].sequence.substr(0,k-1).c_str());
        if (first[0].second) out.emplace_back(first[0].first,nid);
        else in.emplace_back(first[0].first,nid);
        skf_ovl.create_kmers_direction(last,sg.nodes[nid].sequence.substr(sg.nodes[nid].sequence.size()-k+1,k-1).c_str());
        if (last[0].second) in.emplace_back(last[0].first,-nid);
        else out.emplace_back(last[0].first,-nid);
    }
    std::sort(in.begin(),in.end());
    std::sort(out.begin(),out.end());
    //connect out->in for all combinations on each kmer
    uint64_t next_out_idx=0;
    for(auto &i:in){
        while(next_out_idx<out.size() and out[next_out_idx].first<i.first) ++next_out_idx;
        for (auto oidx=next_out_idx;oidx<out.size() and out[oidx].first==i.first;++oidx) {
            sg.add_link(i.second,out[oidx].second,-k+1); //no support, although we could add the DBG operation as such
        }
    }
    sdglib::OutputLog()<<"Graph construction finished"<<std::endl;
}

void GraphMaker::new_graph_from_paired_datastore(const PairedReadsDatastore& ds,  int k, int min_coverage, int num_batches) {
    new_graph_from_kmerlist_128(BatchKmersCounter::countKmersToList(ds, k, min_coverage, num_batches),k);
}

void GraphMaker::new_graph_from_long_datastore(const LongReadsDatastore& ds,  int k, int min_coverage, int num_batches) {
    new_graph_from_kmerlist_128(BatchKmersCounter::countKmersToList(ds, k, min_coverage, num_batches),k);
}