//
// Created by Bernardo Clavijo (EI) on 08/03/2018.
//

#include "FlowFollower.hpp"

Flow FlowFollower::flow_from_node(sgNodeID_t n,float min_winner,float max_looser) {
    auto tags = ws.linked_read_mappers[0].get_node_tags(llabs(n));
    BufferedLRSequenceGetter blrsg(ws.linked_read_datastores[0],1000000,1000);
    auto kmers = ws.linked_read_datastores[0].get_tags_kmers(31, 6, tags,blrsg);
    //std::cout<<"Flowing from "<<n<<" with "<<kmers.size()<<" kmers"<<std::endl;
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

Flow FlowFollower::flow_from_node_kmers(sgNodeID_t n,const std::unordered_set<uint64_t> &kmers) {
    //std::cout<<"Flowing from "<<n<<" with "<<kmers.size()<<" kmers"<<std::endl;
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
    sglib::OutputLog()<<"Create flows starting for "<<nodes.size()<<" nodes."<<std::endl;
    std::vector<sgNodeID_t> nv;
    nv.reserve(nodes.size());
    for (auto &n:nodes)nv.push_back(n);
#pragma omp parallel for schedule(static,1)
    for (auto i=0;i<nv.size();++i){
        auto ff=flow_from_node(nv[i],1,0);
        auto rf=flow_from_node(-nv[i],1,0);
//        if (ff.nodes.size()>2 or rf.nodes.size()>2) {
//            auto ci=ws.kci.compute_compression_for_node(nv[i],1);
//            std::cout<<"Node "<<nv[i]<<" with lenght "<<ws.sg.nodes[nv[i]].sequence.size()<<", "<<
//                     ws.linked_read_mappers[0].get_node_tags(nv[i]).size() << " tags and ci="<<ci
//                     << " produced flows of lenghts "<<ff.nodes.size()<<" and "<<rf.nodes.size()<<std::endl;
//        }
#pragma omp critical
        {
            flows[nv[i]] =ff;
            flows[-nv[i]] =rf;
        }
        //std::cout<<"."<<std::flush;
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
    sglib::OutputLog()<<"DONE!"<<std::endl;
}


void FlowFollower::select_from_all_nodes(uint32_t min_size, uint32_t max_size, uint16_t min_tags, uint16_t max_tags,
                                         float min_ci, float max_ci) {
    uint64_t tnodes=0,ttags=0;
    nodes.clear();
    for (auto n=1;n<ws.sg.nodes.size();++n){
        if (ws.sg.nodes[n].sequence.size()<min_size) continue;
        if (ws.sg.nodes[n].sequence.size()>max_size) continue;
        if (ws.sg.get_fw_links(n).size()!=1 or ws.sg.get_bw_links(n).size()!=1 ) continue;
        auto ntags=ws.linked_read_mappers[0].get_node_tags(n);
        if (ntags.size()<min_tags or ntags.size()>max_tags) continue;
        auto ci=ws.kci.compute_compression_for_node(n,1);
        if (ci<min_ci or ci>max_ci) continue;
        ++tnodes;
        nodes.emplace(n);
        ttags += ntags.size();


    }
    std::cout << "Selected "<<tnodes<<" / "<<ws.sg.nodes.size()<<" with a total "<<ttags<<" tags"<< std::endl;
}

/**
 * @brief sorts nodes so consecutive nodes share tags, then uses a tag kmer cache to speed up flows.
 *
 */
void FlowFollower::create_flows_fast() {
    sglib::OutputLog()<<"Create flows fast starting for "<<nodes.size()<<" nodes into paths collection"<<std::endl;
    //reate an independent buffered tags-kmer-creator per node in a parallel for that creates nodes, saving the
    //flows simply in a vector
    std::vector<Flow> fflows(ws.sg.nodes.size());
    std::vector<Flow> rflows(ws.sg.nodes.size());
    std::atomic<uint64_t> nodesp(0),validflows(0);
    std::vector<sgNodeID_t> nv;
    nv.reserve(nodes.size());
    for (auto &n:nodes)nv.push_back(n);
#pragma omp parallel
    {
        BufferedTagKmerizer btk(ws.linked_read_datastores[0],31,5000,200000,1000);
#pragma omp for schedule(static,100)
        for (auto i=0;i<nv.size();++i){
            auto n=nv[i];
            auto nkmers=btk.get_tags_kmers(6,ws.linked_read_mappers[0].get_node_tags(n));//XXX: hardcoded 2 tags
            fflows[n]=flow_from_node_kmers(n,nkmers);
            rflows[n]=flow_from_node_kmers(-n,nkmers);
//            if (fflows[n].nodes.size()>2 or rflows[n].nodes.size()>2) {
//                auto ci=ws.kci.compute_compression_for_node(n,1);
//                std::cout<<"Node "<<n<<" with lenght "<<ws.sg.nodes[n].sequence.size()<<", "<<
//                         ws.linked_read_mappers[0].get_node_tags(n).size() << " tags and ci="<<ci
//                         << " produced flows of lenghts "<<fflows[n].nodes.size()<<" and "<<rflows[n].nodes.size()<<std::endl;
//            }
            if (fflows[n].nodes.size()>2)++validflows;
            if (rflows[n].nodes.size()>2)++validflows;
            uint64_t np=++nodesp;
            if (np%100==0) sglib::OutputLog()<<np<<" nodes processed, "<<validflows<<" informative flows"<<std::endl;
        }
    }
    sglib::OutputLog()<<" DONE!"<<std::endl;
    //Next create a PathsCollection in the graph
    std::cout<<"Creating paths collection"<<std::endl;
    ws.path_datastores.emplace_back(ws.sg);
    for (auto ff:fflows){
        if (ff.nodes.size()<3) continue;
        ws.path_datastores.back().origin.emplace_back(ff.nodes[0]);
        ws.path_datastores.back().paths.emplace_back(ws.sg,ff.nodes);
    }
    for (auto rf:rflows){
        if (rf.nodes.size()<3) continue;
        ws.path_datastores.back().origin.emplace_back(rf.nodes[0]);
        ws.path_datastores.back().paths.emplace_back(ws.sg,rf.nodes);
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
    std::vector<sgNodeID_t> nv;
    nv.reserve(nodes.size());
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