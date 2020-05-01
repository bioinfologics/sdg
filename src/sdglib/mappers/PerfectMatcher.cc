//
// Created by Bernardo Clavijo (EI) on 30/04/2020.
//

#include "PerfectMatcher.hpp"
#include <sdglib/views/NodeView.hpp>

bool inline are_complement(const char &A,const char &B){
    return (A=='A' and B=='T') or (A=='C' and B=='G') or (A=='G' and B=='C') or (A=='T' and B=='A');
}

void PerfectMatchPart::extend(const std::string & readseq,const std::string & nodeseq) {//this jsut grows until it can't and sets the flags.
    //while(i<size(read) and j<size(node) and read[i]==node[j]) ++i (consider RC, maybe just write the conditions appropriately?)
    auto srp=read_position;
    if (node>0){
        while(read_position<readseq.size()-1 and node_position<nodeseq.size()-1 and nodeseq[node_position+1]==readseq[read_position+1]){
            ++node_position;
            ++read_position;
        }
        completed_node = (node_position==nodeseq.size()-1);
        completed_read = (read_position==readseq.size()-1);
    }
    else {
        while(read_position<readseq.size()-1 and node_position>0 and are_complement(nodeseq[node_position-1],readseq[read_position+1])){
            --node_position;
            ++read_position;
        }
        completed_node = (node_position==0);
        completed_read = (read_position==readseq.size()-1);
    }
    extended=true;
    if (read_position==srp and previous_part!=-1) invalid=true;
}

void PerfectMatchExtender::reset(std::string _readseq){
    matchparts.clear();
    best_path.clear();
    readseq=_readseq;
    last_readpos=0;
}

void PerfectMatchExtender::add_starting_match(sgNodeID_t node_id, uint64_t read_offset, uint64_t node_offset){ //add a new start
    matchparts.emplace_back();
    auto &mp=matchparts.back();
    mp.previous_part=-1;
    mp.node=node_id;
    mp.extended=false;
    mp.read_position=read_offset+k-1;
    mp.extended=mp.completed_read=(mp.read_position==readseq.size()-1);
    if (node_id>0) {
        mp.node_position=node_offset+k-1;
        mp.completed_node = (mp.node_position==dg.sdg.get_node_size(node_id)-1);
    }
    else {
        mp.node_position=node_offset;
        mp.completed_node = (node_offset==0);
    }
}

void PerfectMatchExtender::extend_fw(){
    //TODO: this could be extended backwards in the original hit to add diferentiation on overlap-transitioned hits after an error.
//        std::cout<<"extend_fw called with "<<matchparts.size()<<" starting matchparts"<<std::endl;
    //First, check if any matchparts overlap with each other and start the hit in the OVL, if they do, invalidate the incoming
    for(uint64_t from=0;from<matchparts.size();++from){
        for(uint64_t to=0;to<matchparts.size();++to){
            if (from==to) continue;
            auto d=dg.min_distance(matchparts[from].node,matchparts[to].node);
            if (d<=-k){
                if (matchparts[from].node<0 and matchparts[from].node_position+k<=-d) matchparts[from].invalid=true;
                if (matchparts[from].node>0 and dg.sdg.get_node_size(matchparts[from].node)-matchparts[from].node_position<=-d) matchparts[from].invalid=true;
            }
        }
    }

    for(uint64_t next=0;next<matchparts.size();++next){
        //extend, if end of node add all nexts as unextended parts.
        if (matchparts[next].invalid) continue;
//            std::cout<<"extending matchpart "<<next<<" to node "<<matchparts[next].node<<" with current readpos="<<matchparts[next].read_position<<" and nodepos="<<matchparts[next].node_position<<std::endl;
        matchparts[next].extend(readseq,dg.sdg.nodes[llabs(matchparts[next].node)].sequence);
//            std::cout<<" -> readpos="<<matchparts[next].read_position<<(matchparts[next].completed_read ? " (completed)":"")<<", nodepos="<<matchparts[next].node_position<<(matchparts[next].completed_node ? " (completed)":"")<<std::endl;
        if (matchparts[next].completed_node and not matchparts[next].completed_read){
            for (const auto & l: dg.get_nodeview(matchparts[next].node).next()){
                if (l.distance()>-k+1) continue;
//                    std::cout<<"extending to next node "<<l.node().node_id()<<std::endl;
                matchparts.emplace_back();//add next node part, pointing to this one as previous.
                auto &mp=matchparts.back();
                mp.node=l.node().node_id();
                mp.read_position=matchparts[next].read_position;
                if (mp.node>0) {
                    mp.node_position=-l.distance()-1;
                }
                else {
                    mp.node_position=dg.sdg.get_node_size(mp.node)+l.distance();
                }
                mp.extended=false;
                mp.previous_part=next;
            }
        }
    }
}
void PerfectMatchExtender::set_best_path(){ //set extension_size when returning a path that was extended.
    //find if there is a single part that goes further in the read than all the others.
    uint64_t best_length=0;
    int64_t next_index=0;
    for (auto i=0;i<matchparts.size();++i){
        auto &mp=matchparts[i];
        if (mp.invalid) continue;
        if (mp.read_position==best_length) next_index=-1;
        else if (mp.read_position>best_length) {
            best_length=mp.read_position;
            next_index=i;
        }
    }
    if (next_index==-1){//there was a tie, try to go back to a match that was shared among all paths
        //go through all match-parts, count how many were tied and walk back from those voting.
        //at each stage the first one that has votes from all previous winners is considered the consensus hit so far.
        std::vector<int> votes(matchparts.size());
        int winners_seen=0;
        for (auto i=0;i<matchparts.size();++i) {
            auto &mp = matchparts[i];
            if (mp.invalid or mp.read_position!=best_length) continue;
            ++winners_seen;
            next_index=-1;
            for (auto j=i;j!=-1;j=matchparts[j].previous_part){
                ++votes[j];
                if (votes[j]==winners_seen) next_index=j;
            }
        }
//            if (next_index!=-1) std::cout<<"backtracked to a common part as winner, at index"<<next_index<<std::endl;
    }
    last_nodepos=matchparts[next_index].node_position;


    while (next_index!=-1) {
        best_path.emplace_back(matchparts[next_index].node);
        next_index=matchparts[next_index].previous_part;
    }
    std::reverse(best_path.begin(),best_path.end());
    if (!best_path.empty()) last_readpos=best_length;
    else last_nodepos=0;
}