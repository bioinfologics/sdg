//
// Created by Bernardo Clavijo (EI) on 01/03/2018.
//

#include "TagWalker.hpp"

void TagWalker::remove_crosstalk() {
    for (auto ta:tagsA){
        if (tagsB.count(ta)>0) tags_shared.insert(ta);
    }
    for (auto st:tags_shared) {
        tagsA.erase(st);
        tagsB.erase(st);
    }
    std::cout<<tags_shared.size()<<" tags removed as cross-talk, now A has "<<tagsA.size()<<" and B has "<<tagsB.size()<<" tags"<<std::endl;
}

SequenceGraphPath TagWalker::walk(float min_winner, float max_looser) {
    //option: start from HSPNPs HSPNTs HSPN*, kmerize in any kmer size, find exclussive kmers, get all involved kmers,
    // - Analyse cross-talk
    remove_crosstalk();
    // - Select completely-classified tags, and walk with those.

}