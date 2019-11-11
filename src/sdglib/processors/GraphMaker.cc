//
// Created by Bernardo Clavijo (EI) on 05/07/2018.
//

#include <sdglib/bloom/BloomFilter.hpp>
#include <sdglib/batch_counter/BatchKmersCounter.hpp>
#include "GraphMaker.hpp"

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
void GraphMaker::new_graph_from_kmerset_trivial(const std::unordered_set<uint64_t> & kmerset,uint8_t k) {
    std::cout<<"Constructing Sequence Graph from "<<kmerset.size()<<" "<<std::to_string(k)<<"-mers"<<std::endl;
    std::set<std::pair<uint64_t,uint64_t>> unitig_ends;
    std::vector<std::pair<uint64_t,sgNodeID_t>> kmerovl_bw_nodes,kmerovl_fw_nodes;
    std::unordered_set<uint64_t> used_kmers;
    used_kmers.reserve(kmerset.size());
    //find unitig ends
    for (auto start_kmer:kmerset) {
        if (used_kmers.count(start_kmer)) continue;
        //std::cout<<"processing kmer "<<start_kmer<<": "<<kmer_to_sequence(start_kmer,k)<<std::endl;
        //TODO: can be optimised not to walk the same path many times.
        //start with a kmer, walk fw till neiughbour count !=1
        auto curr_kmer=start_kmer;
        while (true) {
            used_kmers.insert(curr_kmer);
            //first check there is only one kmer going fw
            uint64_t next;
            bool next_found=0;
            for (uint64_t fnk:kmer_fw_neighbours(curr_kmer, k)) {
                if (kmerset.count(kmer_cannonical(fnk,k))) {
                    if (next_found) {
                        next_found=false;
                        break;
                    }
                    else {
                        next=fnk;
                        next_found=true;
                    }
                }
            }
            //then check it only has 1 bw neighbour
            if (next_found) {
                int bc=0;
                for (auto bnk:kmer_bw_neighbours(next,k)) if (kmerset.count(kmer_cannonical(bnk,k))) ++bc;
                if (bc!=1) next_found=false;
            }

            if (!next_found) break;
            if (used_kmers.count(next))break;
            curr_kmer=next;
        }
        auto end_kmer=curr_kmer;
        curr_kmer=start_kmer;
        while (true) {
            used_kmers.insert(curr_kmer);
            //first check there is only one kmer going fw
            uint64_t next;
            bool next_found=0;
            for (uint64_t fnk:kmer_bw_neighbours(curr_kmer, k)) {
                if (kmerset.count(kmer_cannonical(fnk,k))) {
                    if (next_found) {
                        next_found=false;
                        break;
                    }
                    else {
                        next=fnk;
                        next_found=true;
                    }
                }
            }
            //then check it only has 1 bw neighbour
            if (next_found) {
                int bc=0;
                for (auto bnk:kmer_fw_neighbours(next,k)) if (kmerset.count(kmer_cannonical(bnk,k))) ++bc;
                if (bc!=1) next_found=false;
            }

            if (!next_found) break;
            if (used_kmers.count(next))break;
            curr_kmer=next;
        }
        auto begin_kmer=curr_kmer;
        //std::cout<<"begin_kmer: "<<begin_kmer<<std::endl;
        //std::cout<<"end_kmer: "<<end_kmer<<std::endl;
        //add unitig to set
        if (kmer_cannonical(begin_kmer,k)<kmer_cannonical(end_kmer,k)){
            unitig_ends.insert(std::make_pair(begin_kmer,end_kmer));
        }
        else {
            unitig_ends.insert(std::make_pair(kmer_reverse(end_kmer,k),kmer_reverse(begin_kmer,k)));
        }

    }
    std::cout<<"there will be "<<unitig_ends.size()<<" unitigs in the graph"<<std::endl;
    //add each unitig to graph, and its ends to the ovl vector
    char bases[4]={'A','C','G','T'};
    for (auto uends:unitig_ends){
        //create sequence for unitig
        std::string seq=kmer_to_sequence(uends.first,k);
        //add unitig to graph, save id
        auto curr_kmer=uends.first;
        while (curr_kmer!=uends.second){
            auto fwn=kmer_fw_neighbours(curr_kmer,k);
            for (auto i=0;i<4;++i){
                if (kmerset.count(kmer_cannonical(fwn[i],k))){
                    seq.push_back(bases[i]);
                    curr_kmer=fwn[i];
                    break;
                }
            }
        }
        sgNodeID_t nodeid=sg.add_node(Node(seq));
        //std::cout<<"Unitig sequence for node "<<nodeid<<": "<<seq<<std::endl;
        //insert bw overlap of begin kmer (check canonical rep)
        auto begin_ovl=kmer_back_ovl(uends.first,k);
        auto begin_ovl_can=kmer_cannonical(begin_ovl,k-1);
        if (begin_ovl==begin_ovl_can) {
            kmerovl_fw_nodes.push_back(std::make_pair(begin_ovl,nodeid));
            //std::cout<<"FW "<<kmer_to_sequence(begin_ovl,k-1)<<" "<<nodeid<<std::endl;
        }
        else {
            kmerovl_bw_nodes.push_back(std::make_pair(begin_ovl_can,nodeid));
            //std::cout<<"BW "<<kmer_to_sequence(begin_ovl_can,k-1)<<" "<<nodeid<<std::endl;
        }

        //insert fw overlap of end kmer (check canonical rep)
        auto end_ovl=kmer_fw_ovl(uends.second,k);
        auto end_ovl_can=kmer_cannonical(end_ovl,k-1);
        if (end_ovl==end_ovl_can) {
            kmerovl_bw_nodes.push_back(std::make_pair(end_ovl,-nodeid));
            //std::cout<<"BW "<<kmer_to_sequence(end_ovl,k-1)<<" "<<-nodeid<<std::endl;
        }
        else {
            kmerovl_fw_nodes.push_back(std::make_pair(end_ovl_can,-nodeid));
            //std::cout<<"FW "<<kmer_to_sequence(end_ovl_can,k-1)<<" "<<-nodeid<<std::endl;
        }
    }
    //sort ovl vectors
    std::sort(kmerovl_bw_nodes.begin(),kmerovl_bw_nodes.end());
    std::sort(kmerovl_fw_nodes.begin(),kmerovl_fw_nodes.end());
    //for (auto kbn:kmerovl_bw_nodes) std::cout<< "BW "<<kbn.first<<" "<<kbn.second<<std::endl;

    //for (auto kfn:kmerovl_fw_nodes) std::cout<< "FW "<<kfn.first<<" "<<kfn.second<<std::endl;
    //pair nodes in ovl vectors and create links between them (connect from bw to fw when kmers match)
    for (auto kbn:kmerovl_bw_nodes){
        for (auto kfn:kmerovl_fw_nodes){
            if (kbn.first==kfn.first) {
                sg.add_link(kbn.second,kfn.second,-k+1);
            }
        }
    }

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

std::array<__uint128_t,4> kmer_fw_neighbours128_canonical(__uint128_t kmer, uint8_t k) {
    std::array<__uint128_t,4> n;
    for (__uint128_t i=0;i<4;++i){
        n[i]=kmer_cannonical128(((kmer<<2)+i)%(((__uint128_t) 1)<<(k*2)),k);
    }
    return n;
}

std::array<__uint128_t,4> kmer_bw_neighbours128_canonical(__uint128_t kmer, uint8_t k) {
    std::array<__uint128_t,4> n;
    for (__uint128_t i=0;i<4;++i){
        n[i]=kmer_cannonical128((kmer>>2)+(((__uint128_t) i)<<((k-1)*2)),k);
    }
    return n;
}

__uint128_t kmer_back_ovl128(__uint128_t kmer, uint8_t k){
    return kmer>>2;
}

__uint128_t kmer_fw_ovl128(__uint128_t kmer, uint8_t k){
    return  kmer%(((__uint128_t) 1)<<((k-1)*2));
}

void GraphMaker::new_graph_from_kmerset_trivial128(const std::unordered_set<__uint128_t, int128_hash> & kmerset,uint8_t k) {
    //std::cout<<"Constructing Sequence Graph from "<<kmerset.size()<<" "<<std::to_string(k)<<"-mers"<<std::endl;
    std::set<std::pair<__uint128_t,__uint128_t>> unitig_ends;
    std::vector<std::pair<__uint128_t,sgNodeID_t>> kmerovl_bw_nodes,kmerovl_fw_nodes;
    std::unordered_set<__uint128_t, int128_hash> used_kmers;
    used_kmers.reserve(kmerset.size());
    //find unitig ends
    for (__uint128_t start_kmer:kmerset) {
        if (used_kmers.count(start_kmer)) continue;
        __uint128_t curr_kmer=start_kmer;
        while (true) {
            used_kmers.insert(curr_kmer);
            //first check there is only one kmer going fw
            __uint128_t next;
            bool next_found=0;
            for (__uint128_t fnk:kmer_fw_neighbours128(curr_kmer, k)) {
                if (kmerset.count(kmer_cannonical128(fnk,k))) {
                    if (next_found) {
                        next_found=false;
                        break;
                    }
                    else {
                        next=fnk;
                        next_found=true;
                    }
                }
            }
            //then check it only has 1 bw neighbour
            if (next_found) {
                int bc=0;
                for (__uint128_t bnk:kmer_bw_neighbours128(next,k)) if (kmerset.count(kmer_cannonical128(bnk,k))) ++bc;
                if (bc!=1) next_found=false;
            }

            if (!next_found) break;
            if (used_kmers.count(next))break;
            curr_kmer=next;
        }
        auto end_kmer=curr_kmer;
        curr_kmer=start_kmer;
        while (true) {
            used_kmers.insert(curr_kmer);
            //first check there is only one kmer going fw
            __uint128_t next;
            bool next_found=0;
            for (__uint128_t fnk:kmer_bw_neighbours128(curr_kmer, k)) {
                if (kmerset.count(kmer_cannonical128(fnk,k))) {
                    if (next_found) {
                        next_found=false;
                        break;
                    }
                    else {
                        next=fnk;
                        next_found=true;
                    }
                }
            }
            //then check it only has 1 bw neighbour
            if (next_found) {
                int bc=0;
                for (__uint128_t bnk:kmer_fw_neighbours128(next,k)) if (kmerset.count(kmer_cannonical128(bnk,k))) ++bc;
                if (bc!=1) next_found=false;
            }

            if (!next_found) break;
            if (used_kmers.count(next))break;
            curr_kmer=next;
        }
        __uint128_t begin_kmer=curr_kmer;
        //add unitig to set
        if (kmer_cannonical128(begin_kmer,k)<kmer_cannonical128(end_kmer,k)){
            unitig_ends.insert(std::make_pair(begin_kmer,end_kmer));
        }
        else {
            unitig_ends.insert(std::make_pair(kmer_reverse128(end_kmer,k),kmer_reverse128(begin_kmer,k)));
        }

    }
    //std::cout<<"there will be "<<unitig_ends.size()<<" unitigs in the graph"<<std::endl;
    //add each unitig to graph, and its ends to the ovl vector
    char bases[4]={'A','C','G','T'};
    for (auto uends:unitig_ends){
        //create sequence for unitig
        std::string seq=kmer_to_sequence128(uends.first,k);
        auto curr_kmer=uends.first;
        while (curr_kmer!=uends.second){
            auto fwn=kmer_fw_neighbours128(curr_kmer,k);
            for (auto i=0;i<4;++i){
                if (kmerset.count(kmer_cannonical128(fwn[i],k))){
                    seq.push_back(bases[i]);
                    curr_kmer=fwn[i];
                    break;
                }
            }
        }
        //add unitig to graph, save id
        sgNodeID_t nodeid=sg.add_node(Node(seq));
        //insert bw overlap of begin kmer (check canonical rep)
        auto begin_ovl=kmer_back_ovl128(uends.first,k);
        auto begin_ovl_can=kmer_cannonical128(begin_ovl,k-1);
        if (begin_ovl==begin_ovl_can) {
            kmerovl_fw_nodes.push_back(std::make_pair(begin_ovl,nodeid));
        }
        else {
            kmerovl_bw_nodes.push_back(std::make_pair(begin_ovl_can,nodeid));
        }

        //insert fw overlap of end kmer (check canonical rep)
        auto end_ovl=kmer_fw_ovl128(uends.second,k);
        auto end_ovl_can=kmer_cannonical128(end_ovl,k-1);
        if (end_ovl==end_ovl_can) {
            kmerovl_bw_nodes.push_back(std::make_pair(end_ovl,-nodeid));
        }
        else {
            kmerovl_fw_nodes.push_back(std::make_pair(end_ovl_can,-nodeid));
        }
    }
    //sort ovl vectors
    std::sort(kmerovl_bw_nodes.begin(),kmerovl_bw_nodes.end());
    std::sort(kmerovl_fw_nodes.begin(),kmerovl_fw_nodes.end());

    //pair nodes in ovl vectors and create links between them (connect from bw to fw when kmers match)
    for (auto kbn:kmerovl_bw_nodes){
        for (auto kfn:kmerovl_fw_nodes){
            if (kbn.first==kfn.first) {
                sg.add_link(kbn.second,kfn.second,-k+1);
            }
        }
    }

}

class klidxs {
public:
    klidxs(__uint128_t _kmer,uint64_t _idx):kmer(_kmer),idx(_idx){};
    __uint128_t kmer;
    uint64_t idx;
};

inline void get_fw_idxs(std::vector<klidxs> &out, const __uint128_t kmer, uint8_t k,const std::vector<__uint128_t> & kmerlist){
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

inline void get_bw_idxs(std::vector<klidxs> &out, const __uint128_t kmer, uint8_t k,const std::vector<__uint128_t> & kmerlist){
    out.clear();
    __uint128_t cnext;
    //std::vector<__uint128_t>::const_iterator cnitr;
    for (const __uint128_t &bnk:kmer_bw_neighbours128(kmer, k)) {
        cnext=kmer_cannonical128(bnk,k);
        auto cnitr=std::lower_bound(kmerlist.begin(),kmerlist.end(),cnext);
        if (*cnitr==cnext)
            out.emplace_back(bnk,cnitr-kmerlist.begin());
    }
}

bool is_end_bw(const __uint128_t &kmer, uint8_t &k,const std::vector<__uint128_t> & kmerlist){
    std::vector<klidxs> next;
    get_bw_idxs(next,kmer,k,kmerlist);
    if (next.size()!=1) return true;
    auto p=next[0].kmer;
    get_fw_idxs(next,p,k,kmerlist);
    if (next.size()!=1) return true;
    return false;
}

bool is_end_fw(const __uint128_t &kmer, uint8_t &k,const std::vector<__uint128_t> & kmerlist){
    std::vector<klidxs> next;
    get_fw_idxs(next,kmer,k,kmerlist);
    if (next.size()!=1) return true;
    auto p=next[0].kmer;
    get_bw_idxs(next,p,k,kmerlist);
    if (next.size()!=1) return true;
    return false;
}

void GraphMaker::new_graph_from_kmerlist_trivial128(const std::vector<__uint128_t> & kmerlist,uint8_t k) {
    //TODO: add a bloom filter to speed this up
    std::cout<<"Constructing Sequence Graph from "<<kmerlist.size()<<" "<<std::to_string(k)<<"-mers"<<std::endl;
    std::string s;
    s.reserve(1000000); //avoid contig-sequence growth
    std::vector<bool> used_kmers(kmerlist.size());
    std::vector<klidxs> bw,fw;

    for (uint64_t start_kmer_idx=0;start_kmer_idx<kmerlist.size();++start_kmer_idx) {
        if (used_kmers[start_kmer_idx]) continue; //any kmer can only belong to one unitig.ww

        //Check this k-mer is an end/junction
        auto start_kmer = kmerlist[start_kmer_idx];
        bool end_bw = is_end_bw(start_kmer, k, kmerlist);
        auto end_fw = is_end_fw(start_kmer, k, kmerlist);

        if (!end_bw and !end_fw) continue;

        if (end_bw and end_fw) {
            //kmer as unitig
            s=kmer_to_sequence128(start_kmer, k);
            used_kmers[start_kmer_idx] = true;
        } else {
            //unitig start on kmer
            //make sure the unitig starts on FW

            auto current_kmer = start_kmer;
            used_kmers[start_kmer_idx] = true;

            if (end_fw) {
                current_kmer = kmer_reverse128(start_kmer, k);
                end_fw = end_bw;
            }
            //Add kmer as unitig
            s=kmer_to_sequence128(current_kmer, k);

            std::vector<klidxs> fwn;
            unsigned char nucleotides[4] = {'A', 'C', 'G', 'T'};



            //std::cout<<"Starting sequence construction at kmer "<<kmer_to_sequence128(current_kmer,k)<<std::endl;
            while (!end_fw) {
                //Add end nucleotide, update current_kmer;
                get_fw_idxs(fwn, current_kmer, k, kmerlist);
                current_kmer = fwn[0].kmer;
                if (used_kmers[fwn[0].idx]) break; //this should never happen, since these have ends.
                used_kmers[fwn[0].idx]=true;
                s.push_back(nucleotides[current_kmer % 4]);
                end_fw = is_end_fw(current_kmer, k, kmerlist);
            }
        }
        sg.add_node(Node(s));
        if (!sg.nodes.back().is_canonical()) sg.nodes.back().make_rc();//inefficient, but once in a node
    }
    //If there are any perfect circles, they won't have ends, so just pick any unused kmer and go fw till you find the same k-mer fw.
    //(this is going to be even tricker to parallelise)
    std::cout<<"doing the circles now"<<std::endl;
    for (uint64_t start_kmer_idx=0;start_kmer_idx<kmerlist.size();++start_kmer_idx) {
        if (used_kmers[start_kmer_idx]) continue; //any kmer can only belong to one unitig.ww

        //Check this k-mer is an end/junction
        auto start_kmer = kmerlist[start_kmer_idx];
        auto end_fw = is_end_fw(start_kmer, k, kmerlist);
        std::vector<klidxs> fwn;
        unsigned char nucleotides[4] = {'A', 'C', 'G', 'T'};
        auto current_kmer = start_kmer;
        used_kmers[start_kmer_idx] = true;
        s=kmer_to_sequence128(current_kmer, k);

        while (true) {
            //Add end nucleotide, update current_kmer;
            get_fw_idxs(fwn, current_kmer, k, kmerlist);
            current_kmer = fwn[0].kmer;
            if (current_kmer==start_kmer) break;
            used_kmers[fwn[0].idx] = true;
            s.push_back(nucleotides[current_kmer % 4]);
        }
        sg.add_node(Node(s));
        if (!sg.nodes.back().is_canonical()) sg.nodes.back().make_rc();//inefficient, but once in a node
    }

    //save the (k-1)mer in (rev on first k-1 / fw on last k-1) or out ( fw on first k-1 / bw on last k-1)
    std::vector<std::pair<__uint128_t,sgNodeID_t>> in, out;
    in.reserve(2*sg.nodes.size());
    out.reserve(2*sg.nodes.size());
    CStringKMerFactory128 skf_ovl(k-1);
    for (auto nid=1;nid<sg.nodes.size();++nid){
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
}

void GraphMaker::new_graph_from_paired_datastore(const PairedReadsDatastore& ds,  int k, int min_coverage, int num_batches) {
    new_graph_from_kmerlist_trivial128(BatchKmersCounter::countKmersToList(ds, k, min_coverage, num_batches),k);
}

