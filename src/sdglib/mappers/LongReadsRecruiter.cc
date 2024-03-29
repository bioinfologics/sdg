//
// Created by Bernardo Clavijo (EI) on 16/01/2020.
//

#include "sdglib/workspace/WorkSpace.hpp"
#include "LongReadsRecruiter.hpp"
#include "PerfectMatcher.hpp"
#include <sdglib/indexers/NKmerIndex.hpp>
#include <atomic>
#include <sdglib/views/NodeView.hpp>
#include <sdglib/graph/ReadThreadsGraph.hpp>

std::vector<PerfectMatch> PerfectMatchesFilter::truncate_turnaround(const std::vector<PerfectMatch> &in) const {
    if (in.empty()) return {};
    std::unordered_map<sgNodeID_t,uint64_t> node_count;
    std::vector<bool> turn_score;

    //Turn score is true if there's 10% of mathces in opositve direction for this node before this match
    turn_score.reserve(in.size());
    for (auto &m:in){
        ++node_count[m.node];
        turn_score.emplace_back( ((float) node_count[-m.node]) / (node_count[m.node]+node_count[-m.node]) >.1 );
    }

    //Turn fw perc counts how many matches have turn_score true from this match until the end of the read
    std::vector<float> turn_fw_perc;
    turn_fw_perc.reserve(in.size());
    uint64_t f=0,b=0;
    for (auto it=turn_score.crbegin();it<turn_score.crend();++it) {
        if (*it) b+=1;
        else f+=1;
        turn_fw_perc.emplace_back(((float) b)/(b+f));
    }
    std::reverse(turn_fw_perc.begin(),turn_fw_perc.end());

    //If a match is "turned" and >70% of remaining matches are turned too, truncate the read there.
    uint64_t trunc=0;
    while(++trunc<turn_fw_perc.size()){
        if (turn_score[trunc] and turn_fw_perc[trunc]>.7) break;
    }
    return std::vector<PerfectMatch>(in.begin(),in.begin()+trunc);
}

std::vector<PerfectMatch> PerfectMatchesFilter::matches_fw_from_node(sgNodeID_t node, const std::vector<PerfectMatch> &in) const {
    auto last=in.cend();
    for (auto it=in.cbegin();it<in.cend();++it) if (it->node==node) last=it;
    if (last==in.cend()) return {};
    std::vector<PerfectMatch> out(last+1,in.cend());
    if (out.empty()) return {};
    auto offset=last->read_position+ws.sdg.get_node_size(node)-last->node_position;
    for (auto &m:out) m.read_position-=offset;
    return out;
}

std::vector<PerfectMatch> PerfectMatchesFilter::clean_linear_groups(const std::vector<PerfectMatch> &in,int group_size,int small_node_size) const {
    if (in.size()<group_size) return in;
    std::vector<int> mgroups(in.size());//group for each match
    std::vector<int> group_elements(in.size()+1);
    std::vector<int> group_ends(in.size()+1);
    int group=0;

    for (auto i=0;i<in.size();++i){
        if (mgroups[i]) continue;
        mgroups[i]=++group;
        auto group_node=in[i].node;
        auto group_pos=in[i].node_position;
        group_ends[group]=i;
        group_elements[group]=1;
        int hits_since_last=0;
        for (auto j=i+1;j<in.size();++j){
            if (in[j].node!=group_node or in[j].node_position<group_pos){
                if (++hits_since_last==5) break;
            }
            else {
                hits_since_last=0;
                mgroups[j]=group;
                ++group_elements[group];
                group_ends[group]=j;
                group_pos=in[j].node_position;
            }
        }
    }

    int last_big_group=0;
    int last_big_group_end=0;
    for (auto i=0;i<in.size();++i){
        if (i>last_big_group_end) last_big_group=0;
        if (last_big_group==mgroups[i]) continue;
        if (group_elements[mgroups[i]]>=group_size){
            last_big_group=mgroups[i];
            last_big_group_end=group_ends[last_big_group];
        }
        else{
            if (last_big_group) mgroups[i]=0; //reset group for small matches within large groups
        }
    }
    std::vector<PerfectMatch> out;
    out.reserve(in.size());
    for (auto i=0;i<in.size();++i){
        if (mgroups[i]!=0 and (group_elements[mgroups[i]]>=group_size or ws.sdg.get_node_size(in[i].node)<small_node_size))
            out.emplace_back(in[i]);
    }
    return out;
}



std::vector<PerfectMatch> PerfectMatchesFilter::merge_and_sort(const std::vector<std::vector<PerfectMatch> > &in) const {
    //merge first, this assumes matches are actually in distance-from-node form, really.
    //first things first, drop matches to any nodes that do not appear in at least X reads.
    std::unordered_map<sgNodeID_t,uint64_t> node_readcount;
    std::set<sgNodeID_t> nset;
    for (auto &pms:in) {
        nset.clear();
        for (auto &m:pms) nset.insert(m.node);
        for (auto &n:nset) ++node_readcount[n];
    }
    std::vector<std::vector<PerfectMatch> > filtered_in;
    filtered_in.reserve(in.size());
    for (auto &pms:in){
        filtered_in.emplace_back();
        filtered_in.back().reserve(pms.size());
        for (auto &m:pms) if (node_readcount[m.node]>=3) filtered_in.back().push_back(m);
    }

    //TODO: a last bit of cleanup, remove unordered matches to the same node, arising from approximate repeats and sequencing errors


    //now start advancing on each read from the beginning, pick the next read to include using the following criteria:
    // 1 - A read where its next close match to the current node is the one that comes closest next
    // 2 - A read that contains a match to a node that comes before any matches to that node on any other read.
    // 3 - A read that contains the match closest to the current match.

    std::vector<uint64_t> next_in_read(filtered_in.size()); //this keeps tabs on where the reads are
    uint64_t finished_reads=0;
    for (auto i=0;i<filtered_in.size();++i) {
        if (filtered_in[i].empty()){
            next_in_read[i]=-1; //mark
            ++finished_reads;
        }
    }

    std::vector<PerfectMatch> out; //TODO: reserve with the size of all matches?
    sgNodeID_t last_node=0,last_pos=0;
//    uint64_t pass=0;
    int cont=0;
    while(finished_reads<filtered_in.size()){
        //First choose which read to advance
        int64_t winner=-1;
        int64_t winner_index=-1;
        int64_t winner_pos=INT64_MAX;
//        std::cout<<"starting search of the winner"<<std::endl;
        //TODO: periodically check if reads are goind too fast or falling behind and kill them, if too many killed reads, stop

        //TODO: choose a target node to get to -> a node present in most reads (use a distance since last match in each read)

        //TODO: get all the previous nodes out of the way (include them in some order, these will be minoritary nodes anyway)

        //TODO: go through the target node in all reads -> back to start


        while (winner==-1) {
//            if (++pass>200) {
//                std::cout<<"too many passes, aborting..."<<std::endl;
//                break;
//            }
            //First condition: any reads have a match to the last node, plausible of being the next match? FAILS ON LOOPS!!!
            for (auto i = 0; i < filtered_in.size(); ++i) {
                if (next_in_read[i] == -1) continue;
                for (auto j = next_in_read[i]; j < next_in_read[i]+5 and j < filtered_in[i].size(); ++j) {
                    if (filtered_in[i][j].node == last_node) {
                        if (filtered_in[i][j].node_position < winner_pos and filtered_in[i][j].node_position>=last_pos) {
                            //avoid loop problems, by not allowing a massive jump on the read's delta
                            if (next_in_read[i]>0 and filtered_in[i][next_in_read[i]-1].node==last_node and (filtered_in[i][j].node_position-filtered_in[i][next_in_read[i]-1].node_position) * 1.0 / (filtered_in[i][j].read_position-filtered_in[i][next_in_read[i]-1].read_position) < .5 ) {
                                std::cout<<"skipping hit on read "<<i<<" @ "<<j<<" because of big jump!"<<std::endl;
                                continue;
                            }

                            winner = i;
                            winner_index = j;
                            winner_pos = filtered_in[i][j].node_position;
//                                                    std::cout<<"new winner at read "<<i<<" position "<<j<<std::endl;
                        }
                        //break;
                    }
                }
            }
            //If no reads can advance on the last node, try to guess waht the most immediate node is and set it as last_node
            if (winner==-1){
                bool ln_set=false;
                //Pick a node that no other read have as preceeded by other nodes.
                std::cout<<"trying to find a front-only node "<<std::endl;
                std::unordered_map<sgNodeID_t,uint64_t> candidates;
                for (auto i = 0; i < filtered_in.size(); ++i) {
                    if (next_in_read[i]==-1) continue;
                    ++candidates[filtered_in[i][next_in_read[i]].node];
                }


                //A node can be in 3 "status": front, absent, later, ideally we want a node that is only front and absent
                //create a set of all nodes that front, then go through every read
                std::unordered_map<sgNodeID_t,uint64_t> nexts;
                for (auto i = 0; i < filtered_in.size(); ++i) {
                    if (next_in_read[i] == -1) continue;
                    std::unordered_set<sgNodeID_t> next_nodes;
                    auto firstnode=filtered_in[i][next_in_read[i]].node;
                    for (auto j = next_in_read[i]+1; j < next_in_read[i] + 10 and j < filtered_in[i].size(); ++j){
                        if (filtered_in[i][j].node!=firstnode) next_nodes.insert(filtered_in[i][j].node);
                    }
                    for (auto &nn:next_nodes) ++nexts[nn];
                }
                std::cout<<"candidates are:";
                for (auto &c:candidates)std::cout<<"    "<<c.first<<" ( "<<c.second<<" - "<<nexts[c.first]<<" )";
                std::cout<<std::endl;

                //first choice: the most-voted node with no nexts attached
                auto max_first=0;
                for (auto &c:candidates) {
                    if (nexts[c.first]==0 and c.second>max_first) {
                        max_first=c.first;
                        last_node=c.first;
                        last_pos=0;
                        ln_set=true;
                    }
                }
                if (not ln_set) {

                    //second choice: the node with the least nexts
                    auto min_nexts = filtered_in.size();
                    for (auto &c:candidates) {
                        if (nexts[c.first] < min_nexts) {
                            min_nexts = nexts[c.first];
                            last_node = c.first;
                            last_pos = 0;
                        }
                    }
                }
                std::cout<<"could not find a continuation from last_node, trying with node "<<last_node<<std::endl;
                if (++cont>3) break;
            }
            else cont=0;
        }
//        if (++pass>1000) break;
        //Now advance the read, set last_node and last_pos
        //std::cout<<"read "<<winner<<" wins with the match at position "<<winner_index<<std::endl;
        auto &m=filtered_in[winner][winner_index];
        std::cout<<"WINNER read #"<<winner<<" @ "<<winner_index<<" -> "<<m.node<<" @ "<<m.node_position<<std::endl;
        for (;next_in_read[winner]<=winner_index;++next_in_read[winner]) {
            auto &m=filtered_in[winner][next_in_read[winner]];
            std::cout<<"inserting read #"<<winner<<" @ "<<next_in_read[winner]<<" -> "<<m.node<<" @ "<<m.node_position<<std::endl;
            out.emplace_back(filtered_in[winner][next_in_read[winner]]);
        }
//        std::cout<<"after matches insertion there is "<<out.size()<<" matches in the output collection"<<std::endl;
        last_node=out.back().node;
        last_pos=out.back().node_position;

        if (next_in_read[winner]==filtered_in[winner].size()) {
            next_in_read[winner]=-1;
            ++finished_reads;
        }
    }
    return out;
}

void PerfectMatchesMergeSorter::init_from_node(sgNodeID_t n, const LongReadsRecruiter & lrr, int min_reads, int group_size, int small_node_size) {
    out.clear();

    PerfectMatchesFilter pmf(ws);
    //get all reads in the node
    read_matches.clear();
    read_matches.reserve(lrr.node_reads[llabs(n)].size());
    std::unordered_map<sgNodeID_t,uint64_t> node_read_count;
    std::unordered_set<sgNodeID_t> read_nodes;
    //add matches, reversed if need be, cleaned and fw from node
    for (auto rid:lrr.node_reads[llabs(n)]){
        if (n<0) rid=-rid;
        read_matches.emplace_back( pmf.matches_fw_from_node(n, pmf.truncate_turnaround( pmf.clean_linear_groups(
            rid>0 ? lrr.read_perfect_matches[rid] :
            lrr.reverse_perfect_matches(lrr.read_perfect_matches[-rid],lrr.datastore.read_to_fileRecord[-rid].record_size)
            ,group_size,small_node_size))));

        if (read_matches.back().empty()) read_matches.pop_back();
        else {
            read_nodes.clear();
            for (auto &m:read_matches.back()) read_nodes.insert(m.node);
            for (auto &n:read_nodes) ++node_read_count[n];
        }
    }
    //remove any nodes not in enough reads
    for (auto &rm:read_matches){
        rm.erase(std::remove_if(rm.begin(),rm.end(),[&node_read_count,&min_reads](PerfectMatch &m){return node_read_count[m.node]<min_reads;}),rm.end());
    }
    //remove any newly empty reads due to node filtering
    read_matches.erase(std::remove_if(read_matches.begin(), read_matches.end(),[](std::vector<PerfectMatch> &v){return v.empty();}),read_matches.end());

    read_last_hit_position.resize(read_matches.size(),0);
    read_dropped_position.resize(read_matches.size(),-1);
    read_next_match.resize(read_matches.size(),-0);

}

void PerfectMatchesMergeSorter::find_next_node(int d, float candidate_percentaje, float first_percentaje, bool verbose) {
    next_node=0;
    //explore the next x bp of reads, mark nodes appearing there.
    std::unordered_set<sgNodeID_t> read_nodes;
    std::unordered_map<sgNodeID_t,uint64_t> node_read_count;
    //count in how many reads each node appears, only use next hits up to last_hit_position + d
    uint64_t used_reads=0;
    for (auto i=0;i<read_matches.size();++i) {
        if (read_next_match[i]==-1) continue;
        ++used_reads;
        read_nodes.clear();
        for (auto it=read_matches[i].cbegin()+read_next_match[i];
            it<read_matches[i].cend() and it->read_position-read_last_hit_position[i]<d;
            ++it)
            read_nodes.insert(it->node);
        for (auto &n:read_nodes) ++node_read_count[n];
    }
    //any node that appears in 80% of the reads is a safe node to get to, check which one of them comes first and check any nodes that come before
    //TODO: consider that some reads just won't be long enough to get to the node! Compute distance to first match of node in read and only use reads that get there
    if (verbose) std::cout<<used_reads<<" reads have hits to candidate nodes"<<std::endl;
    std::unordered_set<sgNodeID_t> popular_nodes;
    for (auto &nc:node_read_count){
        if (nc.second>used_reads*candidate_percentaje) {
            popular_nodes.insert(nc.first);
        }
    }
    //pick the earliest of the popular nodes
    std::unordered_map<sgNodeID_t,uint64_t> first_node_read_count;
    uint64_t total_reads_with_first=0;
    for (auto i=0;i<read_matches.size();++i) {
        if (read_next_match[i] == -1) continue;
        for (auto it=read_matches[i].cbegin()+read_next_match[i];
             it<read_matches[i].cend() and it->read_position-read_last_hit_position[i]<d;
             ++it)
            if (popular_nodes.count(it->node)) {
                ++first_node_read_count[it->node];
                ++total_reads_with_first;
                break;
            }
    }
    for (auto fnc:first_node_read_count){
        if (verbose) std::cout<<"Popular node "<<fnc.first<<" is the first of the popular nodes in "<<fnc.second<<" / "<< node_read_count[fnc.first] <<" reads "<<std::endl;
        if (fnc.second>first_percentaje*node_read_count[fnc.first]) {
            if (next_node!=0) {
                if (verbose) std::cout<<"there's more than one potential next node, this version of the algorithm can't deal with that"<<std::endl;
                next_node=0;
                break;
            }
            next_node=fnc.first;
        }
    }

}

void PerfectMatchesMergeSorter::advance_reads_to_node() {
    //TODO: as of now this only skips the intermediate nodes, sorry!
    /*std::unordered_map<sgNodeID_t,uint64_t> pre_nodes;
    std::unordered_set<sgNodeID_t> seen_pre_nodes;
    for (auto i=0;i<read_matches.size();++i) {
        if (read_next_match[i] == -1) continue;
        seen_pre_nodes.clear();
        for (auto it = read_matches[i].cbegin() + read_next_match[i];
             it < read_matches[i].cend();
             ++it) {
            //std::cout << "Read " << i << " hits node " << it->node << " before getting to next node" << std::endl;
            if (it->node == next_node) {
                for (auto &pn:seen_pre_nodes) ++pre_nodes[pn];
                break;
            }
            seen_pre_nodes.insert(it->node);
        }
    }
    for (auto &pn:pre_nodes){
        //std::cout<<"Pre node "<<pn.first<<" is present in "<<pn.second<<" reads"<<std::endl;
    }*/
    //TODO: compute a reasonable distance to the node and drop all other reads
    for (auto i = 0; i < read_matches.size(); ++i) {
        if (read_next_match[i] == -1) continue;
        for (auto it = read_matches[i].cbegin() + read_next_match[i];
             it < read_matches[i].cend();
             ++it) {
            //std::cout << "Read " << i << " hits node " << it->node << " before getting to next node" << std::endl;
            if (it->node == next_node) {
                read_next_match[i]=it-read_matches[i].cbegin();
                break;
            }
        }
    }
}
void PerfectMatchesMergeSorter::advance_reads_through_node() {

    int64_t last_pos=0;

    while (true) {
        int winner=-1;
        int64_t winner_index=0;
        int64_t winner_pos=INT64_MAX;
        for (auto i = 0; i < read_matches.size(); ++i) {
            if (read_next_match[i] == -1) continue;
            for (auto j = read_next_match[i]; j < read_next_match[i] + 5 and j < read_matches[i].size(); ++j) {
                if (read_matches[i][j].node == next_node) {
                    if (read_matches[i][j].node_position < winner_pos and
                        read_matches[i][j].node_position >= last_pos) {
                        //avoid loop problems, by not allowing a massive jump on the read's delta
                        if (read_next_match[i] > 0 and read_matches[i][read_next_match[i] - 1].node == next_node and
                            (read_matches[i][j].node_position - read_matches[i][read_next_match[i] - 1].node_position) *
                            1.0 /
                            (read_matches[i][j].read_position - read_matches[i][read_next_match[i] - 1].read_position) <
                            .5) {
                            //std::cout << "skipping hit on read " << i << " @ " << j << " because of big jump!" << std::endl;
                            continue;
                        }

                        winner = i;
                        winner_index = j;
                        winner_pos = read_matches[i][j].node_position;
                        //std::cout<<"new winner at read "<<i<<" position "<<j<<std::endl;
                    }
                    //break;
                }
            }
        }
        if (winner==-1) break;

        auto &m=read_matches[winner][winner_index];
        //std::cout<<"WINNER read #"<<winner<<" @ "<<winner_index<<" -> "<<m.node<<" @ "<<m.node_position<<std::endl;
        for (;read_next_match[winner]<=winner_index;++read_next_match[winner]) {
            auto &m=read_matches[winner][read_next_match[winner]];
            //std::cout<<"inserting read #"<<winner<<" @ "<<read_next_match[winner]<<" -> "<<m.node<<" @ "<<m.node_position<<std::endl;
            out.emplace_back(read_matches[winner][read_next_match[winner]]);
        }
//        std::cout<<"after matches insertion there is "<<out.size()<<" matches in the output collection"<<std::endl;
        last_pos=out.back().node_position;
        read_last_hit_position[winner]=out.back().read_position;

        if (read_next_match[winner]==read_matches[winner].size()) {
            read_next_match[winner]=-1;
            //++finished_reads;
        }

    }
}

void PerfectMatchesMergeSorter::drop_conflictive_reads() {

}

LongReadsRecruiter::LongReadsRecruiter(SequenceDistanceGraph &_sdg, const LongReadsDatastore &datastore,uint8_t _k, uint16_t _f):
    sdg(_sdg),datastore(datastore),k(_k),f(_f) {
    reset_recruitment();
};


void LongReadsRecruiter::dump(std::string filename) {
    std::ofstream ofs(filename);
    sdglib::write_flat_vectorvector(ofs,node_reads);
    sdglib::write_flat_vectorvector(ofs,read_perfect_matches);
}

void LongReadsRecruiter::dump_threads(std::string filename) {
    std::ofstream ofs(filename);
    sdglib::write_flat_vectorvector(ofs,read_threads);
}

void LongReadsRecruiter::load(std::string filename) {
    std::ifstream ifs(filename);
    sdglib::read_flat_vectorvector(ifs,node_reads);
    sdglib::read_flat_vectorvector(ifs,read_perfect_matches);
}

void LongReadsRecruiter::load_threads(std::string filename) {
    std::ifstream ifs(filename);
    sdglib::read_flat_vectorvector(ifs,read_threads);
}

void LongReadsRecruiter::reset_recruitment() {
    node_reads.clear();
    node_reads.resize(sdg.nodes.size());
    read_perfect_matches.clear();
    read_perfect_matches.resize(datastore.size()+1);
}
//TODO: optimise the in-vector class to do a binary search (right now the kmatches are not well-ordered because of RCs)
template<class T>
inline bool in_vector(const std::vector<T> &V,T VAL) { return std::find(V.begin(),V.end(),VAL)!=V.end();}

void LongReadsRecruiter::perfect_mappings(uint16_t seed_size, uint64_t first_read, uint64_t last_read) {
    if (seed_size<k) seed_size=k;
    NKmerIndex nki(sdg,k,f);
    if (last_read==0) last_read=datastore.size();
#pragma omp parallel shared(nki)
    {
        StreamKmerFactory skf(k);
        std::vector<std::pair<bool,uint64_t >> read_kmers; //TODO: type should be defined in the kmerisers
        std::vector<std::vector<std::pair<int32_t,int32_t >>> kmer_matches; //TODO: type should be defined in the index class
        PerfectMatch pmatch(0,0,0);
        ReadSequenceBuffer sequence_reader(datastore);
#pragma omp for
        for (auto rid=first_read;rid<=last_read;++rid){
            std::vector<PerfectMatch> private_read_perfect_matches;
            auto seq=sequence_reader.get_read_sequence(rid);
            read_kmers.clear();
            skf.produce_all_kmers(seq,read_kmers);
            nki.get_all_kmer_matches(kmer_matches,read_kmers);
            std::map<sgNodeID_t,uint32_t > node_match_count;
            auto rksize=read_kmers.size();
            auto longest_seed = seed_size;
            for (auto i=0;i<rksize;++i) {
                longest_seed=(longest_seed>=seed_size ? longest_seed-1:longest_seed); //decreases the length of the skipped matches in each cycle
                pmatch.node = 0;
                for (auto const &m:kmer_matches[i]) {
                    auto last_node = m.first;
                    auto last_pos = m.second;
                    auto needed_hit=i+longest_seed-k;
                    if (needed_hit>i and (needed_hit >= rksize or not in_vector(kmer_matches[needed_hit], {last_node, last_pos + needed_hit -i})))
                        continue;
                    auto last_i = i;
                    while (last_i + 1 < rksize and
                           in_vector(kmer_matches[last_i + 1], {last_node, last_pos + last_i - i + 1}))
                        ++last_i;
                    if (last_i == needed_hit and i==pmatch.read_position) {
                        pmatch.node = 0;
                    } else if (last_i > needed_hit) {
                        pmatch.node = last_node;
                        pmatch.node_position = (last_node>0 ? last_pos:last_pos+1-k);
                        pmatch.read_position = i;
                        longest_seed = pmatch.size = k + last_i - i;
                    }
                }
                if (pmatch.node != 0) {
                    private_read_perfect_matches.emplace_back(pmatch);
                }
            }
            read_perfect_matches[rid] = private_read_perfect_matches;
        }
    }
}

std::vector<PerfectMatch> LongReadsRecruiter::reverse_perfect_matches(const std::vector<PerfectMatch> &matches,uint64_t rsize) const {
    if (matches.empty()) return {};
    for (const auto &m:matches) if (rsize<m.read_position+m.size) rsize=m.read_position+m.size;

    std::vector<PerfectMatch> rmatches;
    rmatches.reserve(matches.size());

    for (int i=matches.size()-1;i>=0;--i){
        const auto &m=matches[i];
        rmatches.emplace_back(-m.node,sdg.get_node_size(m.node)-m.node_position-m.size,rsize-m.read_position-m.size,m.size);
    }
    std::sort(rmatches.begin(),rmatches.end());
    return rmatches;
}

void LongReadsRecruiter::map(uint16_t seed_size) {
    if (seed_size<k) seed_size=k;
    NKmerIndex nki(sdg,k,f);
    sdglib::OutputLog()<<"Index created with "<<nki.end()-nki.begin()<<" items"<<std::endl;
    map_to_index(seed_size,nki);
}

void LongReadsRecruiter::anchormap(uint16_t seed_size, std::string kcname, std::string countname, uint32_t fmin,
                                   uint32_t fmax, uint32_t graph_fmin, uint32_t graph_fmax) {
    std::unordered_set<uint64_t> whitelist;
    auto &kc=sdg.ws.get_kmer_counter(kcname);
    if (kc.k!=k) throw std::runtime_error("LRR anchormap: k is different between KC and LRR");
    auto kcg=kc.get_count_by_name("sdg");
    auto kcr=kc.get_count_by_name(countname);
    sdglib::OutputLog()<<"Creating whitelist from "<<kc.kindex.size()<<" kmers"<<std::endl;
    for (auto kcix=0;kcix<kc.kindex.size();++kcix){
        if (kcg[kcix]>=graph_fmin and kcg[kcix]<=graph_fmax and kcr[kcix]>=fmin and kcr[kcix]<=fmax)
            whitelist.insert(kc.kindex[kcix]);
    }
    if (seed_size<k) seed_size=k;
    sdglib::OutputLog()<<"Creating index with "<<whitelist.size()<<" kmers in whitelist"<<std::endl;
    NKmerIndex nki(sdg,k,whitelist);
    sdglib::OutputLog()<<"Index created with "<<nki.end()-nki.begin()<<" items"<<std::endl;

    map_to_index(seed_size,nki);
}

void LongReadsRecruiter::map_to_index(uint16_t seed_size, NKmerIndex &nki) {
    auto total_reads = datastore.size();
    auto maped_reads_count=0;
    bool debug=false;
#pragma omp parallel shared(nki)
    {
        StreamKmerFactory skf(k);
        std::vector<std::pair<int64_t,uint64_t >> read_kmers; //TODO: type should be defined in the kmerisers
        std::vector<std::vector<std::pair<int32_t,int32_t >>> kmer_matches; //TODO: type should be defined in the index class
        PerfectMatchExtender pme(sdg,k);
        ReadSequenceBuffer sequence_reader(datastore);
        uint64_t private_mapped_reads=0;
#pragma omp for schedule(dynamic,1)
        for (auto rid=1;rid<=datastore.size();++rid){

            std::vector<PerfectMatch> private_read_perfect_matches;
            const auto seq=std::string(sequence_reader.get_read_sequence(rid));
            if (debug) std::cout<<"Read "<< rid << " --> " << seq << std::endl;
            read_kmers.clear();
            skf.produce_all_kmers(seq.c_str(),read_kmers);
            pme.set_read(seq);
            int64_t last_end=0;
            for (const auto &rpk:read_kmers){
                //coordinates in the rpk start at 1/-1, and are basically the number of the k-mer in the sequence, so
                if (llabs(rpk.first)+k<=last_end+2) {
                    if (debug) std::cout<<"skipping k-mer at "<<rpk.first<<" because llabs("<<rpk.first<<")-1+"<<(int)k<<"<"<<last_end<<std::endl;
                    continue; //skip until the k-mer doesn't belong to last hit
                }
                if (debug) std::cout<<"Next unmatched kmer at rp "<< rpk.first<<std::endl;
                auto kmatch=nki.find(rpk.second);
                if (kmatch!=nki.end() and kmatch->kmer==rpk.second){

                    pme.reset();
                    for (;kmatch!=nki.end() and kmatch->kmer==rpk.second;++kmatch) {
                        if (debug) std::cout<<"Match found! to "<<kmatch->contigID<<"@"<<kmatch->offset<<std::endl;
                        auto contig=kmatch->contigID;
                        int64_t pos=kmatch->offset-1;
                        if (pos<0) {
                            contig=-contig;
                            pos=-pos-2;
                        }
                        if (rpk.first<0) contig=-contig;
                        pme.add_starting_match(contig, llabs(rpk.first)-1, pos);
                    }
                    pme.extend_fw();
                    pme.set_best_path();
                    const auto & pmebp=pme.best_path;
                    if (debug) std::cout<<"After extension last_readpos is "<<pme.last_readpos<<std::endl;
                    //last_readpos is 0-based and left to the end USED base, while rpk is 1-based to allow for sign as direction
                    //hence, an alignment at the very beginning of the read would be rpk.first=1 and last_readpos=30
                    if (pme.last_readpos-llabs(rpk.first)>=seed_size-2) {
                        if (!pmebp.empty()) {
                            last_end=pme.last_readpos; //avoid extra index lookups for kmers already used once
                            if (debug) std::cout<<"Matches being used, so last_end moved to "<<last_end<<std::endl;
                            pme.make_path_as_perfect_matches();
                            private_read_perfect_matches.insert(private_read_perfect_matches.end(),pme.best_path_matches.begin(),pme.best_path_matches.end());
                            if (debug) std::cout<<pme.best_path_matches.size()<<" matches from best path kept"<<std::endl;
                        }
                    }
//                    //TODO: matches shold be extended left to avoid unneeded indeterminaiton when an error occurrs in an overlap region and the new hit matches a further part of the genome.
//                    std::cout<<"rki after match extension: "<<rki<<" / "<<read_kmers.size()<<std::endl;
                }

            }
            read_perfect_matches[rid] = private_read_perfect_matches;
            if (++private_mapped_reads%5000==0) {
                maped_reads_count+=5000;
                sdglib::OutputLog()<<private_mapped_reads<<" reads mapped on thread "<<omp_get_thread_num() << ". Progress: " << maped_reads_count << "/"<< total_reads <<std::endl;
            }
        }
    }
}

void LongReadsRecruiter::clean_sandwich_matches(int look_ahead) {
    std::vector<PerfectMatch> cleaned;
    for (auto &rpm:read_perfect_matches){
        if (rpm.size()<3) continue;
        cleaned.clear();
        cleaned.push_back(rpm[0]);
        for (auto i=1;i<rpm.size()-1;++i) {

            if (rpm[i].node==cleaned.back().node) cleaned.push_back(rpm[i]);
            else {
                bool skip=false;
                for (auto j=i+1;j<rpm.size() and j<=i+look_ahead;++j) {
                    auto &nm=rpm[j];
                    if (nm.node==rpm[i].node) break;
                    if (nm.node==cleaned.back().node and nm.node_position>cleaned.back().node_position) {
                        skip=true;
                        break;
                    }
                }
                if (not skip) cleaned.push_back(rpm[i]);
            }


        }
        cleaned.push_back(rpm.back());
        rpm=cleaned;
    }
}

//TODO: add read position and node position to match
void LongReadsRecruiter::recruit_reads(uint16_t seed_size, uint16_t seed_count, int64_t first_read,
                                       int64_t last_read) {
    node_reads.clear();
    node_reads.resize(sdg.nodes.size());
    if (last_read==0) last_read=datastore.size();
    for (int64_t rid=first_read;rid<=last_read;++rid) {
        std::map<sgNodeID_t, uint32_t> node_match_count;
        for (const auto &pmatch:read_perfect_matches[rid]) {
            if (pmatch.size>=seed_size) {
                node_match_count[pmatch.node] = node_match_count[pmatch.node] + 1;
                if (node_match_count[pmatch.node] == seed_count) {
                    node_reads[llabs(pmatch.node)].emplace_back(pmatch.node>0?rid:-rid);
                }
            }
        }
    }
}
void LongReadsRecruiter::recruit_threads() {
    node_threads.clear();
    node_threads.resize(sdg.nodes.size());
    for (auto rid=1;rid<=read_threads.size();++rid)
        for (const auto &match:read_threads[rid])
                    node_threads[llabs(match.node)].emplace_back(rid);
}

std::vector<std::pair<int,int>> LongReadsRecruiter::find_all_valid_blocks(const std::vector<PerfectMatch> & matches, int min_count, int min_bp){
    std::vector<std::pair<int,int>> blocks;
    std::unordered_map<sgNodeID_t,int> node_count,node_bp;
    for (auto &m:matches) {
        ++node_count[m.node];
        node_bp[m.node]+=m.size;
    }
    int start=0;
    int nid=0;
    int count;
    for (auto i=0;i<matches.size();++i) {
        //std::cout<<matches[i].node<<" "<<node_count[matches[i].node]<<std::endl;
        if (node_count[matches[i].node]<min_count or node_bp[matches[i].node]<min_bp) continue;
        if (matches[i].node!=nid) {
            nid=matches[i].node;
            start=i;
            count=0;
        }
        ++count;
        //std::cout<<nid<<" reached count "<<count<<std::endl;
        if (node_count[nid]==count) blocks.emplace_back(start,i);
    }
    return blocks;
}

void LongReadsRecruiter::simple_thread_reads(int min_count,int min_bp) {
    read_threads.clear();
    read_threads.resize(read_perfect_matches.size());
    for (auto rid=1;rid<read_perfect_matches.size();++rid){
        int64_t thread_offset=0;//offset between thread and read, computed on last block's position vs. its last match
        const auto & matches=read_perfect_matches[rid];
        for (auto nb:find_all_valid_blocks(matches,min_count,min_bp)) {
            int64_t pstart=matches[nb.first].read_position-matches[nb.first].node_position+thread_offset;
            int64_t pend=pstart+sdg.get_node_size(matches[nb.first].node);
            thread_offset=pstart+matches[nb.second].node_position-matches[nb.second].read_position;
            read_threads[rid].emplace_back(matches[nb.first].node, pstart, pend);
        }
    }
}

ReadThreadsGraph LongReadsRecruiter::rtg_from_threads(bool remove_duplicated, int min_thread_nodes) {
    ReadThreadsGraph rtg(sdg);
    std::vector<uint64_t> node_counts(sdg.nodes.size());
    for (auto const &t:read_threads) {
        for (auto const &np:t ) ++node_counts[llabs(np.node)];
    }
    rtg.links.resize(node_counts.size());
    for (uint64_t i=1;i<node_counts.size();++i) rtg.links[i].reserve(2*node_counts[i]);
    for (auto rid=1;rid<read_threads.size();++rid) {
        rtg.add_thread(rid,read_threads[rid],remove_duplicated,min_thread_nodes);
    }
    return rtg;
}

std::vector<sgNodeID_t> LongReadsRecruiter::path_fw(seqID_t read_id, sgNodeID_t node) const {
    //check for obvious no-path conditions:
    if (read_paths[llabs(read_id)].size()<2) return {};
    //get read paths in the right orientation - leaves r2p empty if not using pairs
    std::vector<sgNodeID_t> r1p;
    if (read_id > 0) {
        r1p=read_paths[read_id];
    }
    else {
        sdglib::reverse_path(r1p,read_paths[-read_id]); //no read pair, as the pair comes "BEFORE"
    }

    //discards up to node in r1p
    auto r1it=r1p.cbegin();
    while (r1it<r1p.cend() and *r1it!=node) ++r1it;
    if (r1it==r1p.cend()) return {};
    std::vector<sgNodeID_t> path_fw;
    path_fw.insert(path_fw.end(),++r1it,r1p.cend());

    return path_fw;
}

std::vector<std::vector<sgNodeID_t> > LongReadsRecruiter::all_paths_fw(sgNodeID_t node) const {
    std::vector<std::vector<sgNodeID_t> > r;
    std::unordered_set<seqID_t> rids;
    for (auto rid:node_paths[llabs(node)] ){
        if (node<0) rid=-rid;
        r.emplace_back(path_fw(rid,node));
        if (r.back().empty()) r.pop_back();

    }
    return r;
}
