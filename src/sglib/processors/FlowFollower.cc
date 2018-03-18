//
// Created by Bernardo Clavijo (EI) on 08/03/2018.
//

#include "FlowFollower.hpp"

Flow FlowFollower::flow_from_node(sgNodeID_t n,float min_winner,float max_looser) {
    //std::cout << "Starting walk on " << n << "... " << std::flush;
    auto tags = ws.linked_read_mappers[0].get_node_tags(llabs(n));
    //std::cout << tags.size() << " tags... " << std::flush;
    BufferedLRSequenceGetter blrsg(ws.linked_read_datastores[0],1000000,1000);
    auto kmers = ws.linked_read_datastores[0].get_tags_kmers(31, 6, tags,blrsg);
    //std::cout << kmers.size() << " kmers." << std::endl;
    SequenceGraphPath p(ws.sg, {n});
//    for (auto i = 0; i < 2; ++i) {
        while (true) {
            auto fwl = ws.sg.get_fw_links(p.nodes.back());
            if (fwl.empty()) break;
            std::vector<sgNodeID_t> fw_nodes;
            bool no_distinctive=false;
            for (auto l:fwl) fw_nodes.push_back(l.dest);
            auto distinctivekmers = get_distinctive_kmers(fw_nodes);
            for (auto &dk:distinctivekmers) {
                if (dk.size()==0){
                    //std::cout<<"WARNING: no distinctive kmers for at least one node: ";
                    //for (auto i=0;i<fw_nodes.size();++i) std::cout<<fw_nodes[i]<<"("<<distinctivekmers[i].size()<<") ";
                    //std::cout<<std::endl;
                    no_distinctive=true;
                }
            }
            if (no_distinctive) break;
            sgNodeID_t best = 0, second = 0;
            uint64_t best_score = 100000, second_score = 1000000;

            for (auto i = 0; i < fw_nodes.size(); ++i) {
                std::unordered_set<uint64_t> inters;
                auto &lkmers = distinctivekmers[i];
                uint64_t hits=0,misses = 0;
                for (auto x:lkmers) {
                    if (kmers.count(x)) ++hits;
                    else ++misses;
                }
                if (misses<best_score) {
                    second = best;
                    best = fw_nodes[i];
                    second_score = best_score;
                    best_score = misses;
                } else if (misses < second_score) {
                    second = fw_nodes[i];
                    second_score = misses;
                }

            }
            if (best_score != 0 or second_score==0) {
                //std::cout << "stopping because of score uncertainty " << best_score << " - " << second_score << std::endl;
                break;
            }
            bool b = false;
            for (auto x:p.nodes) if (llabs(x) == llabs(best)) b = true;
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
    std::vector<std::unordered_set<uint64_t>> distinctive_kmers,temp_distinctive;

    if (nodes.size()==1) {
        distinctive_kmers.emplace_back();
        for (auto x:node_kmers[0]) distinctive_kmers.back().insert(x);
    }
    else {
        for (auto nk:node_kmers) {
            distinctive_kmers.emplace_back();
            temp_distinctive.emplace_back();
            for (auto x:nk) {
                if (kmer_counts[x] == 1) {
                    temp_distinctive.back().insert(x);

                } else if (kmer_counts[x] == nodes.size()) {
                    for (auto y:temp_distinctive.back()) distinctive_kmers.back().insert(y);
                    temp_distinctive.back().clear();
                }
            }
            if (distinctive_kmers.back().size()==0 and temp_distinctive.back().size()>0) recover_short=true;
        }
        if (recover_short) {
            uint64_t mind=temp_distinctive[0].size();
            for (auto td:temp_distinctive) if (td.size()<mind) mind=td.size();
            distinctive_kmers.clear();
            for (auto nk:node_kmers) {
                distinctive_kmers.emplace_back();
                for (auto x:nk) {
                    if (kmer_counts[x] == 1) {
                        distinctive_kmers.back().insert(x);
                        if (distinctive_kmers.back().size()>=mind) break;
                    }
                }
            }
        }

    }

    return distinctive_kmers;
}

void FlowFollower::create_flows() {
    std::vector<sgNodeID_t> nv(nodes.size());
    for (auto &n:nodes)nv.push_back(n);
#pragma omp parallel for schedule(static,1)
    for (auto i=0;i<nv.size();++i){
        auto ff=flow_from_node(nv[i],1,0);
        auto rf=flow_from_node(-nv[i],1,0);
#pragma omp critical
        {
            flows[nv[i]] =ff;
            flows[-nv[i]] =rf;
        }
        std::cout<<"."<<std::flush;
        /**std::cout<<"flows["<<n<<"]= ";
        for (auto m:flows[n].nodes) std::cout<<m<<", ";
        std::cout<<std::endl;
        std::cout<<" PATH: ";
        for (auto m:flows[n].nodes) std::cout<<"seq"<<llabs(m)<<", ";
        std::cout<<std::endl;
        std::cout<<"flows["<<-n<<"]= ";
        for (auto m:flows[-n].nodes) std::cout<<m<<", ";
        std::cout<<std::endl;
        std::cout<<" PATH: ";
        for (auto m:flows[-n].nodes) std::cout<<"seq"<<llabs(m)<<", ";
        std::cout<<std::endl;**/
    }
    std::cout<<std::endl;
}

/**
 * @brief sorts nodes so consecutive nodes share tags, then uses a tag kmer cache to speed up flows.
 *
 */
void FlowFollower::create_flows_all_fast() {
    struct node_tag_t{
        sgNodeID_t node;
        //uint16_t unused_rounds=0;
        std::set<bsg10xTag> tags;
        float shared_perc(const struct node_tag_t & other){
            uint64_t shared=0;
            auto i1=tags.begin();
            auto i2=other.tags.begin();
            while (i1!=tags.end() and i2!=other.tags.end()){
                if (*i1==*i2) {
                    ++shared;
                    ++i1;
                    ++i2;
                }
                else if (*i1<*i2) ++i1;
                else ++i2;
            }
            return ((float) shared) / tags.size();
        };
        uint64_t shared(const std::set<bsg10xTag> &td ){
            uint64_t shared=0;
            for (auto t1=tags.begin(),t2=td.begin();t1!=tags.end() and t2!=td.end();){
                if (*t1==*t2) {
                    ++shared;
                    ++t1;
                    ++t2;
                }
                else if (*t1<*t2) ++t1;
                else ++t2;
            }
            return shared;
        };
        uint64_t shared_est(const std::set<bsg10xTag> &td ){
            uint64_t shared=0;
            auto t1lim=tags.size()/10;
            size_t t1c=0;
            for (auto t1=tags.begin(),t2=td.begin();t1c<t1lim and t2!=td.end();){
                if (*t1==*t2) {
                    ++shared;
                    ++t1;
                    ++t1c;
                    ++t2;
                }
                else if (*t1<*t2) {
                    ++t1;
                    ++t1c;
                }
                else ++t2;
            }
            return 10*shared;
        };
    };
    uint64_t tnodes=0,ttags=0;
    std::vector<struct node_tag_t> node_tags;
    node_tags.reserve(ws.sg.nodes.size());
    for (auto n=1;n<ws.sg.nodes.size();++n){
        auto ntags=ws.linked_read_mappers[0].get_node_tags(n);
        if (ntags.size()>=50) {
            ++tnodes;
            node_tags.emplace_back();
            node_tags.back().node=n;
            for (auto &t:ntags) node_tags.back().tags.insert(t);
            ttags+=node_tags.back().tags.size();
        }
    }
    std::cout << "There are "<<tnodes<<" / "<<ws.sg.nodes.size()<<" nodes with 50+ tags, totalling "<<ttags<<" tags"<< std::endl;


    for (auto perc:{.25}) {
        std::vector<bsg10xTag> cached_tags;
        std::set<bsg10xTag > cached_tags_set;
        for (auto &t:node_tags[0].tags) {
            cached_tags.push_back(t);
            cached_tags_set.insert(t);
        }
        uint64_t cache_size=2000;
        uint64_t total_uses=node_tags[0].tags.size(),total_reads=node_tags[0].tags.size(),total_swaps=0;
        uint64_t last_offset=1;
        for (uint64_t i = 1; i < node_tags.size(); ++i) {
            //std::cout<<"Finding good candidates for position "<<i<<std::endl;
            uint64_t j;
            for (auto jb = 0; jb < node_tags.size()-i; ++jb) {
                j = last_offset+jb;
                //std::cout<<"jb="<<jb<<" j="<<j<<" node_tags.size()="<<node_tags.size()<<std::endl;
                if (j>=node_tags.size()) {
                    //std::cout<<"j = i + j - node_tags.size() = "<<i<<" + "<<j<<" - "<<node_tags.size()<<" = ";
                    j=i+j-node_tags.size();
                    //std::cout<<j<<std::endl;
                }
                if (j<i) {
                    std::cout<<"j<i, how the ??"<<std::endl;
                    std::cout<<"i="<<i<<" last_offset="<<last_offset<<" jb="<<jb<<" node_tags.size()="<<node_tags.size()<<" j="<<j<<std::endl;
                }
                auto shared_count = node_tags[j].shared_est(cached_tags_set);
                if (shared_count >= perc * node_tags[j].tags.size()) {
                    //std::cout<<"Found close node at "<<j<<" ("<<shared_count<<"/"<<node_tags[j].tags.size()<<")"<<std::endl;
                    if (j != i) {
                        ++total_swaps;
                        std::swap(node_tags[i].node, node_tags[j].node);
                        std::swap(node_tags[i].tags, node_tags[j].tags);
                    }
                    break;
                }

            }
            last_offset = (j>i? j: i+1);
            std::vector<bsg10xTag> new_cache;
            new_cache.reserve(cache_size);
            new_cache.insert(new_cache.end(), node_tags[i].tags.begin(), node_tags[i].tags.end());
            for (auto &t:cached_tags) {
                if (new_cache.size() >= cache_size) break;
                if (std::find(node_tags[i].tags.begin(), node_tags[i].tags.end(), t) == node_tags[i].tags.end()) {
                    new_cache.push_back(t);
                }
            }
            auto shared_count = node_tags[i].shared(cached_tags_set);
            std::swap(cached_tags, new_cache);
            cached_tags_set.clear();
            for (auto &t:cached_tags) cached_tags_set.insert(t);

            total_uses += node_tags[i].tags.size();
            total_reads += node_tags[i].tags.size() - shared_count;

            if (0 == i % 1000) std::cout << "Position " << i << " : " << total_reads << " / " << total_uses << std::endl;
        }
        std::cout << "After " << total_swaps << " swaps, there will be " << total_reads << " reads for a total of "
                  << total_uses << " uses" << std::endl;
    }

}

/**
 * @brief starts from node, goes + first, - later. Incorporates all transversed flows while deciding. Turns back.
 * @return
 */
SequenceGraphPath FlowFollower::skate_from_node(sgNodeID_t n) {
    struct used_flow_t{Flow flow; bool active; uint32_t pos;} ;
    std::vector<struct used_flow_t> used_flows;
    SequenceGraphPath path(ws.sg,{n});

    std::cout<<"Starting to skate from node "<<n<<std::endl;
    for (auto pass=0;pass<2;++pass) {
        while (true) {
            auto fwls = path.get_next_links();
            if (fwls.empty()) break;
            std::cout<<"Adding node "<<path.nodes.back()<<std::endl;
            if (flows.count(path.nodes.back())>0 and flows[path.nodes.back()].nodes.size()>3) {
                used_flows.push_back({flows[path.nodes.back()], true, 0});
                std::cout<<"Node has a flow with "<<flows[path.nodes.back()].nodes.size()<<" nodes, added as #"<<used_flows.size()-1<<std::endl;
            }
            sgNodeID_t next = 0;
            for (auto i=0;i<used_flows.size();++i) {
                auto &f=used_flows[i];
                //check if flow active
                if (!f.active) continue;
                if (f.pos == f.flow.nodes.size() - 1) {
                    f.active = false;
                    std::cout<<"Flow #"<<i<<" closed"<<std::endl;
                    continue;
                }

                ++f.pos;
                //, option==0 or option==this flow, check if flow is finished
                //raise conflict if flow different (and stop)
                if (0 == next) {
                    next = f.flow.nodes[f.pos];
                } else if (next != f.flow.nodes[f.pos]) {
                    next = 0;
                    break;
                }

            }
            if (next == 0) break;
            bool b=false;
            for (auto x:path.nodes) if (llabs(x)==llabs(next)) b=true;
            if (b) break;
            path.nodes.push_back(next);
        }
        path.reverse();
        used_flows.clear();
    }
    return path;
}

std::vector<SequenceGraphPath> FlowFollower::skate_from_all(int min_node_flow, uint64_t min_path_length) {
    std::vector<SequenceGraphPath> r;
    std::vector<sgNodeID_t> nv(nodes.size());
    for (auto &n:nodes)nv.push_back(n);
//#pragma omp parallel for schedule(static,1)
    for (auto i=0;i<nv.size();++i){
        auto n=nv[i];
        if (flows[n].nodes.size()>=min_node_flow or flows[-n].nodes.size()>=min_node_flow){
            auto p=skate_from_node(n);
            if (p.get_sequence().size()>min_path_length) {
                if (llabs(p.nodes.front())>llabs(p.nodes.back())) p.reverse();
//#pragma omp critical
                {
                    r.push_back(p);

                    std::cout << "PATH skated from " << n << std::endl;
                    for (auto x:p.nodes) std::cout << "seq" << llabs(x) << ",";
                    std::cout << std::endl;
                }
            }
        }
    }
    std::cout<<"Sorting paths."<<std::endl;
    std::sort(r.begin(),r.end());
    r.erase(std::unique(r.begin(),r.end()),r.end());
    std::cout<<"Skate from all finished."<<std::endl;
    return r;
}