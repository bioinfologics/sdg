//
// Created by Bernardo Clavijo (EI) on 08/03/2018.
//

#include "FlowFollower.hpp"

Flow FlowFollower::flow_from_node(sgNodeID_t n,float min_winner,float max_looser) {


    //std::cout << "Starting walk on " << n << "... " << std::flush;
    auto tags = ws.linked_read_mappers[0].get_node_tags(llabs(n));
    //std::cout << tags.size() << " tags... " << std::flush;
    BufferedLRSequenceGetter blrsg(ws.linked_read_datastores[0],1000000,1000);
    auto kmers = ws.linked_read_datastores[0].get_tags_kmers(31, 3, tags,blrsg);
    //std::cout << kmers.size() << " kmers." << std::endl;
    SequenceGraphPath p(ws.sg, {n});
//    for (auto i = 0; i < 2; ++i) {
        while (true) {
            auto fwl = ws.sg.get_fw_links(p.nodes.back());
            if (fwl.empty()) break;
            std::vector<sgNodeID_t> fw_nodes;
            for (auto l:fwl) fw_nodes.push_back(l.dest);
            auto distinctivekmers = get_distinctive_kmers_truncated(fw_nodes);//XXX: now using truncated it should be a first pass with normal, a second with truncated if needed
            for (auto &dk:distinctivekmers) {
                if (dk.size()==0){
                    std::cout<<"WARNING: no distinctive kmers for at least one node: ";
                    for (auto i=0;i<fw_nodes.size();++i) std::cout<<fw_nodes[i]<<"("<<distinctivekmers[i].size()<<") ";
                    std::cout<<std::endl;
                }
            }

            sgNodeID_t best = 0, second = 0;
            double best_score = 0, second_score = 0;

            for (auto i = 0; i < fw_nodes.size(); ++i) {
                std::unordered_set<uint64_t> inters;
                auto &lkmers = distinctivekmers[i];
                double score;
                uint64_t hits = 0;
                for (auto x:lkmers) if (kmers.count(x)) ++hits;
                score = (double) hits / lkmers.size();
                //std::cout << "scoring transition to " << fw_nodes[i] << ": " << hits << "/" << lkmers.size() << "="
                //          << score << std::endl;
                if (best == 0 or score > best_score) {
                    second = best;
                    best = fw_nodes[i];
                    second_score = best_score;
                    best_score = score;
                } else if (second == 0 or score > second_score) {
                    second = fw_nodes[i];
                    second_score = score;
                }
            }
            //std::cout << best_score << " - " << second_score << std::endl;
            if (best_score == 0 or best_score < min_winner or second_score > max_looser) {
                std::cout << "stopping because of score uncertainty "<< best_score << " - " << second_score << std::endl;
                break;
            }
            bool b = false;
            for (auto n:p.nodes) if (llabs(n) == llabs(best)) b = true;
            if (b) {
                //std::cout << "stopping on circular path" << std::endl;
                break;
            }
            p.nodes.push_back(best);
        }
//        p.reverse();
//    }
    Flow f;
    f.nodes=p.nodes;
    return f;
}


std::vector<std::unordered_set<uint64_t>> FlowFollower::get_distinctive_kmers(std::vector<sgNodeID_t> nodes) {
    std::unordered_set<uint64_t> seen_kmers,shared_kmers;
    std::vector<std::unordered_set<uint64_t>> distinctive_kmers;
    for (auto n:nodes){
        distinctive_kmers.emplace_back();
        StringKMerFactory skf(ws.sg.nodes[llabs(n)].sequence,31);
        std::vector<uint64_t> nkmers;
        nkmers.reserve(ws.sg.nodes[llabs(n)].sequence.size());
        skf.create_kmers(nkmers);
        for (auto x:nkmers) {
            if (seen_kmers.count(x) > 0) {
                shared_kmers.insert(x);
            } else {
                distinctive_kmers.back().insert(x);
                seen_kmers.insert(x);
            }
        }
    }
    for (auto &dk:distinctive_kmers) for (auto sk:shared_kmers) if (dk.count(sk)) dk.erase(sk);
    return distinctive_kmers;
}

std::string DNArc(const std::string & s ){
    std::string rc;
    rc.reserve(s.size());
    for (auto b=s.rbegin();b!=s.rend();++b) {
        if ('A'==*b) rc.push_back('T');
        else if ('C'==*b) rc.push_back('G');
        else if ('G'==*b) rc.push_back('C');
        else if ('T'==*b) rc.push_back('A');
        else rc.push_back('N');
    }
    return rc;
}
/**
 * @brief does the same as get_distinctive_kmers, but truncates the nodes' sequences on the last shared k-mer to avoid size-related effects
 * @param nodes
 * @return
 */
std::vector<std::unordered_set<uint64_t>> FlowFollower::get_distinctive_kmers_truncated(std::vector<sgNodeID_t> nodes) {
    //TODO: consider the case where both things are completely different!
    std::vector<std::vector<uint64_t>> node_kmers;
    std::map<uint64_t,uint16_t> kmer_counts;
    for (auto n:nodes){
        node_kmers.emplace_back();
        auto s=(n>0 ? ws.sg.nodes[llabs(n)].sequence:DNArc(ws.sg.nodes[llabs(n)].sequence));
        StringKMerFactory skf(s,31);
        node_kmers.back().reserve(ws.sg.nodes[llabs(n)].sequence.size());
        skf.create_kmers(node_kmers.back());
        for (auto x:node_kmers.back()) {
            if (kmer_counts.count(x)==0) kmer_counts[x]=1;
            else ++kmer_counts[x];
        }
    }
    bool recover_short=false;
    //now go through the nodes:
    //  if the kmers are not in the "shared" set add to temp list
    //  if kmers in "shared_by_all" list, move temp list to distinctive collection.
    std::vector<std::unordered_set<uint64_t>> distinctive_kmers;
    if (nodes.size()==1) {
        distinctive_kmers.emplace_back();
        for (auto x:node_kmers[0]) distinctive_kmers.back().insert(x);
    }
    else {
        for (auto nk:node_kmers) {
            std::vector<uint64_t> temp_distinctive;
            distinctive_kmers.emplace_back();
            for (auto x:nk) {
                if (kmer_counts[x] == 1) {
                    temp_distinctive.push_back(x);

                } else if (kmer_counts[x] == nodes.size()) {
                    for (auto y:temp_distinctive) distinctive_kmers.back().insert(y);
                    temp_distinctive.clear();
                }
            }
            //if (distinctive_kmers.back().size()==0 and temp_distinctive.size()>0) recover_short=true;
        }
    }
    //if (recover_short) {
    //
    //}

    return distinctive_kmers;
}

void FlowFollower::create_flows() {
    for (auto n:nodes){
        flows[n]=flow_from_node(n,1,0);
        std::cout<<"flows["<<n<<"]= ";
        for (auto m:flows[n].nodes) std::cout<<m<<", ";
        std::cout<<std::endl;
        std::cout<<" PATH: ";
        for (auto m:flows[n].nodes) std::cout<<"seq"<<llabs(m)<<", ";
        std::cout<<std::endl;
        flows[-n]=flow_from_node(-n,1,0);
        std::cout<<"flows["<<-n<<"]= ";
        for (auto m:flows[-n].nodes) std::cout<<m<<", ";
        std::cout<<std::endl;
        std::cout<<" PATH: ";
        for (auto m:flows[-n].nodes) std::cout<<"seq"<<llabs(m)<<", ";
        std::cout<<std::endl;

    }
}