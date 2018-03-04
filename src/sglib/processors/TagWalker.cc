//
// Created by Bernardo Clavijo (EI) on 01/03/2018.
//

#include "TagWalker.hpp"

float TagWalker::remove_crosstalk() {
    for (auto ta:tagsA){
        if (tagsB.count(ta)>0) tags_shared.insert(ta);
    }
    for (auto st:tags_shared) {
        tagsA.erase(st);
        tagsB.erase(st);
    }
    //std::cout<<tags_shared.size()<<" tags removed as cross-talk, now A has "<<tagsA.size()<<" and B has "<<tagsB.size()<<" tags"<<std::endl;
    return (tags_shared.size()>0 ? ((float) tags_shared.size()) / (tags_shared.size()+tagsA.size()+tagsB.size()) : 0);
}

void TagWalker::dump_reads(std::string prefix) {
    std::ofstream Ar1of(prefix+"_A_R1.fasta");
    std::ofstream Ar2of(prefix+"_A_R2.fasta");
    std::ofstream Br1of(prefix+"_B_R1.fasta");
    std::ofstream Br2of(prefix+"_B_R2.fasta");
    BufferedLRSequenceGetter blrsg(ws.linked_read_datastores[0],10000000,1000);
    for (auto tag:tagsA) {
        for (auto r:ws.linked_read_datastores[0].get_tag_reads(tag)){
            *(r%2==1 ? &Ar1of:&Ar2of)<<">"<<tag<<"_"<<(r-1)/2<<std::endl<<blrsg.get_read_sequence(r)<<std::endl;
        }
    }
    for (auto tag:tagsB) {
        for (auto r:ws.linked_read_datastores[0].get_tag_reads(tag)){
            *(r%2==1 ? &Br1of:&Br2of)<<">"<<tag<<"_"<<(r-1)/2<<std::endl<<blrsg.get_read_sequence(r)<<std::endl;
        }
    }
}

std::vector<std::unordered_set<uint64_t>> TagWalker::get_distinctive_kmers(std::vector<sgNodeID_t> nodes) {
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

std::vector<SequenceGraphPath> TagWalker::walk(float min_winner, float max_looser) {
    //option: start from HSPNPs HSPNTs HSPN*, kmerize in any kmer size, find exclussive kmers, get all involved kmers,
    for(auto pass=0;pass<2;++pass) {

        sgNodeID_t n= (pass==0?nodeA:nodeB);

        //std::cout << "Starting walk on " << n << "... " << std::flush;
        auto tags = ws.linked_read_mappers[0].get_node_tags(llabs(n));
        //std::cout << tags.size() << " tags... " << std::flush;
        BufferedLRSequenceGetter blrsg(ws.linked_read_datastores[0],1000000,1000);
        auto kmers = ws.linked_read_datastores[0].get_tags_kmers(31, 3, tags,blrsg);
        //std::cout << kmers.size() << " kmers." << std::endl;
        SequenceGraphPath p(ws.sg, {n});
//        {
//            auto lkmers = get_distinctive_kmers({n})[0];
//            double score;
//            uint64_t hits = 0;
//            for (auto x:lkmers) if (kmers.count(x)) ++hits;
//            score = (double) hits / lkmers.size();
//            //std::cout << "Self-score: " << hits << "/" << lkmers.size() << "=" << score << std::endl;
//        }
        for (auto i = 0; i < 2; ++i) {
            while (true) {
                auto fwl = ws.sg.get_fw_links(p.nodes.back());
                if (fwl.empty()) break;
                std::vector<sgNodeID_t> fw_nodes;
                for (auto l:fwl) fw_nodes.push_back(l.dest);
                auto distinctivekmers = get_distinctive_kmers(fw_nodes);


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
                    //std::cout << "stopping because of score uncertainty" << std::endl;
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
            p.reverse();
            (pass==0?pathA:pathB).nodes=p.nodes;
        }
    }
    return {pathA,pathB};
}