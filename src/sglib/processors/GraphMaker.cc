//
// Created by Bernardo Clavijo (EI) on 05/07/2018.
//

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
void GraphMaker::new_graph_from_kmerset_trivial(std::unordered_set<uint64_t> kmerset,uint8_t k) {
    std::cout<<"Constructing Sequence Graph from "<<kmerset.size()<<" "<<std::to_string(k)<<"-mers"<<std::endl;
    std::set<std::pair<uint64_t,uint64_t>> unitig_ends;
    std::vector<std::pair<uint64_t,sgNodeID_t>> kmerovl_bw_nodes,kmerovl_fw_nodes;

    //find unitig ends
    for (auto start_kmer:kmerset) {
        //std::cout<<"processing kmer "<<start_kmer<<": "<<kmer_to_sequence(start_kmer,k)<<std::endl;
        //TODO: can be optimised not to walk the same path many times.
        //start with a kmer, walk fw till neiughbour count !=1
        auto curr_kmer=start_kmer;
        while (true) {
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
            curr_kmer=next;
        }
        auto end_kmer=curr_kmer;
        curr_kmer=start_kmer;
        while (true) {
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
        std::cout<<"Unitig sequence for node "<<nodeid<<": "<<seq<<std::endl;
        //insert bw overlap of begin kmer (check canonical rep)
        auto begin_ovl=kmer_back_ovl(uends.first,k);
        auto begin_ovl_can=kmer_cannonical(begin_ovl,k-1);
        if (begin_ovl==begin_ovl_can) {
            kmerovl_fw_nodes.push_back(std::make_pair(begin_ovl,nodeid));
            std::cout<<"FW "<<kmer_to_sequence(begin_ovl,k-1)<<" "<<nodeid<<std::endl;
        }
        else {
            kmerovl_bw_nodes.push_back(std::make_pair(begin_ovl_can,nodeid));
            std::cout<<"BW "<<kmer_to_sequence(begin_ovl_can,k-1)<<" "<<nodeid<<std::endl;
        }

        //insert fw overlap of end kmer (check canonical rep)
        auto end_ovl=kmer_fw_ovl(uends.second,k);
        auto end_ovl_can=kmer_cannonical(end_ovl,k-1);
        if (end_ovl==end_ovl_can) {
            kmerovl_bw_nodes.push_back(std::make_pair(end_ovl,-nodeid));
            std::cout<<"BW "<<kmer_to_sequence(end_ovl,k-1)<<" "<<-nodeid<<std::endl;
        }
        else {
            kmerovl_fw_nodes.push_back(std::make_pair(end_ovl_can,-nodeid));
            std::cout<<"FW "<<kmer_to_sequence(end_ovl_can,k-1)<<" "<<-nodeid<<std::endl;
        }
    }
    std::cout<<"-----"<<std::endl;
    //sort ovl vectors
    std::sort(kmerovl_bw_nodes.begin(),kmerovl_bw_nodes.end());
    std::sort(kmerovl_fw_nodes.begin(),kmerovl_fw_nodes.end());
    for (auto kbn:kmerovl_bw_nodes) std::cout<< "BW "<<kbn.first<<" "<<kbn.second<<std::endl;

    for (auto kfn:kmerovl_fw_nodes) std::cout<< "FW "<<kfn.first<<" "<<kfn.second<<std::endl;
    //pair nodes in ovl vectors and create links between them (connect from bw to fw when kmers match)
    for (auto kbn:kmerovl_bw_nodes){
        for (auto kfn:kmerovl_fw_nodes){
            if (kbn.first==kfn.first) {
                sg.add_link(kbn.second,kfn.second,-k+1);
            }
        }
    }

}