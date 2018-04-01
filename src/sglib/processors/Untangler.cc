//
// Created by Bernardo Clavijo (EI) on 26/02/2018.
//

#include "Untangler.hpp"
#include "TagWalker.hpp"

struct Counter
{
    struct value_type { template<typename T> value_type(const T&) { } };
    void push_back(const value_type&) { ++count; }
    size_t count = 0;
};

template<typename T1, typename T2>
size_t intersection_size(const T1& s1, const T2& s2)
{
    Counter c;
    set_intersection(s1.begin(), s1.end(), s2.begin(), s2.end(), std::back_inserter(c));
    return c.count;
}

uint64_t Untangler::solve_canonical_repeats_by_tags(std::unordered_set<uint64_t> &reads_to_remap) {
    std::unordered_set<sgNodeID_t> unsolved_repeats;
    unsolved_repeats.clear();
    uint64_t aa_count=0,ab_count=0,unsolved_count=0,non_evaluated=0;

    std::cout << " Finding trivial repeats to analyse with tags" << std::endl;
    std::vector<bool> used(ws.sg.nodes.size());
    std::vector<SequenceGraphPath> paths_solved;
    //if (ws.verbose_log!="") verbose_log_file<<"==== Round of repetition analysis started ===="<<std::endl;
    for (auto n = 1; n < ws.sg.nodes.size(); ++n) {
        if (used[n]) {
            //if (verbose_log!="") verbose_log_file<<n<<" is used"<<std::endl;
            continue;
        }
        auto fwl = ws.sg.get_fw_links(n);
        auto bwl = ws.sg.get_bw_links(n);
        if (fwl.size() != 2 or bwl.size() != 2) {
            //if (verbose_log!="") verbose_log_file<<n<<" has "<<bwl.size()<<" ins and "<<fwl.size()<<" outs"<<std::endl;
            continue;
        }
        auto f0 = fwl[0].dest;
        auto f1 = fwl[1].dest;
        auto b0 = bwl[0].dest;
        auto b1 = bwl[1].dest;
        if (used[(b0 > 0 ? b0 : -b0)] or used[(b1 > 0 ? b1 : -b1)] or used[(f0 > 0 ? f0 : -f0)] or
            used[(f1 > 0 ? f1 : -f1)]) {
            non_evaluated++;
            //if (verbose_log!="") verbose_log_file<<n<<" neighbors are used"<<std::endl;
            continue;
        }
        sgNodeID_t all[4] = {f0, f1, b0, b1};
        bool ok = true;
        for (auto x:all) if (x == n or x == -n or ws.sg.nodes[(x > 0 ? x : -x)].sequence.size() < 199) ok = false;
        for (auto j = 0; j < 3; ++j)
            for (auto i = j + 1; i < 4; ++i)
                if (all[i] == all[j] or all[i] == -all[j])ok = false; //looping node
        if (!ok) {
            //if (verbose_log!="") verbose_log_file<<n<<" and neighbors are short"<<std::endl;
            continue;
        }

//                std::cout<<"Repeat: [ "<<-b0<<" | "<<-b1<<" ] <-> "<<n<<" <-> [ "<<f0<<" | "<<f1<<" ]"<<std::endl;
        std::set<prm10xTag_t> b0tags, b1tags, f0tags, f1tags;

//                std::cout<<"Reads in nodes -> f0:"<<scaff.lrmappers[0].reads_in_node[(f0 > 0 ? f0 : -f0)].size()
//                        <<"   f1:"<<scaff.lrmappers[0].reads_in_node[(f1 > 0 ? f1 : -f1)].size()
//                        <<"   b0:"<<scaff.lrmappers[0].reads_in_node[(b0 > 0 ? b0 : -b0)].size()
//                        <<"   b1:"<<scaff.lrmappers[0].reads_in_node[(b1 > 0 ? b1 : -b1)].size()<<std::endl;

        for (auto rm:ws.linked_read_mappers[0].reads_in_node[(f0 > 0 ? f0 : -f0)])
            f0tags.insert(ws.linked_read_mappers[0].datastore.get_read_tag(rm.read_id));
        for (auto rm:ws.linked_read_mappers[0].reads_in_node[(f1 > 0 ? f1 : -f1)])
            f1tags.insert(ws.linked_read_mappers[0].datastore.get_read_tag(rm.read_id));
        for (auto rm:ws.linked_read_mappers[0].reads_in_node[(b0 > 0 ? b0 : -b0)])
            b0tags.insert(ws.linked_read_mappers[0].datastore.get_read_tag(rm.read_id));
        for (auto rm:ws.linked_read_mappers[0].reads_in_node[(b1 > 0 ? b1 : -b1)])
            b1tags.insert(ws.linked_read_mappers[0].datastore.get_read_tag(rm.read_id));

        f0tags.erase(0);
        f1tags.erase(0);
        b0tags.erase(0);
        b1tags.erase(0);

        std::set<prm10xTag_t> shared1,shared2;

        for (auto t:f0tags) if (f1tags.count(t)>0) {shared1.insert(t);};
        for (auto t:shared1){f0tags.erase(t);f1tags.erase(t);};
        for (auto t:b0tags) if (b1tags.count(t)>0) {shared2.insert(t);};
        for (auto t:shared2) {b0tags.erase(t);b1tags.erase(t);};

        std::set<prm10xTag_t> aa, bb, ba, ab;
        std::set_intersection(b0tags.begin(), b0tags.end(), f0tags.begin(), f0tags.end(),
                              std::inserter(aa, aa.end()));
        std::set_intersection(b0tags.begin(), b0tags.end(), f1tags.begin(), f1tags.end(),
                              std::inserter(ab, ab.end()));
        std::set_intersection(b1tags.begin(), b1tags.end(), f0tags.begin(), f0tags.end(),
                              std::inserter(ba, ba.end()));
        std::set_intersection(b1tags.begin(), b1tags.end(), f1tags.begin(), f1tags.end(),
                              std::inserter(bb, bb.end()));
//                std::cout<<"Tags in   f0: "<<f0tags.size()<<"  f1: "<<f1tags.size()<<"  b0: "<<b0tags.size()<<"  b1: "<<b1tags.size()<<std::endl;
//                std::cout<<"Tags support   aa: "<<aa.size()<<"  bb: "<<bb.size()<<"  ab: "<<ab.size()<<"  ba: "<<ba.size();
        if (aa.size() > 3 and bb.size() > 3 and
            std::min(aa.size(), bb.size()) > 10 * std::max(ab.size(), ba.size())) {
            //std::cout << " Solved as AA BB !!!" << std::endl;
            used[(b0 > 0 ? b0 : -b0)] = true;
            used[(b1 > 0 ? b1 : -b1)] = true;
            used[n] = true;
            used[(f0 > 0 ? f0 : -f0)] = true;
            used[(f1 > 0 ? f1 : -f1)] = true;
            paths_solved.push_back(SequenceGraphPath(ws.sg, {-b0, n, f0}));
            paths_solved.push_back(SequenceGraphPath(ws.sg, {-b1, n, f1}));
            ++aa_count;

        } else if (ba.size() > 3 and ab.size() > 3 and
                   std::min(ba.size(), ab.size()) > 10 * std::max(aa.size(), bb.size())) {
            //std::cout << " Solved as AB BA !!!" << std::endl;
            used[(b0 > 0 ? b0 : -b0)] = true;
            used[(b1 > 0 ? b1 : -b1)] = true;
            used[n] = true;
            used[(f0 > 0 ? f0 : -f0)] = true;
            used[(f1 > 0 ? f1 : -f1)] = true;
            paths_solved.push_back(SequenceGraphPath(ws.sg, {-b0, n, f1}));
            paths_solved.push_back(SequenceGraphPath(ws.sg, {-b1, n, f0}));
            ++ab_count;
        }
        else {
            unsolved_repeats.insert(n);
            ++unsolved_count;
            //if (verbose_log!="") verbose_log_file<<n<<" unsolved aa: "<<aa.size()<<"  bb: "<<bb.size()<<"  ab: "<<ab.size()<<"  ba: "<<ba.size()<<std::endl;
        }
    }

    std::cout << paths_solved.size() << " paths to join" << std::endl;
    //if (paths_solved.size() > 0) mod = true;
    reads_to_remap;
    for (auto &p:paths_solved) {
        ws.sg.join_path(p);
        for (auto n:p.nodes) {
            for (auto rm:ws.linked_read_mappers[0].reads_in_node[(n > 0 ? n : -n)]) {
                reads_to_remap.insert((rm.read_id % 2 ? rm.read_id : rm.read_id - 1));
                ws.linked_read_mappers[0].read_to_node[rm.read_id] = 0;
            }
        }

    }
    std::cout<<"Path analysis summary AA:"<<aa_count<<" AB:"<<ab_count
             <<" Unsolved:"<<unsolved_count<<" Non evaluated:"<<non_evaluated<<std::endl;
    return paths_solved.size();
}

uint64_t Untangler::expand_canonical_repeats_by_tags(float min_ci, float max_ci) {
    uint64_t aa_count=0,ab_count=0,unsolved_count=0;

    std::cout << " Finding trivial repeats to analyse with tags" << std::endl;
    for (auto n = 1; n < ws.sg.nodes.size(); ++n) {
        auto fwl = ws.sg.get_fw_links(n);
        auto bwl = ws.sg.get_bw_links(n);
        if (fwl.size() != 2 or bwl.size() != 2) continue;
        auto f0 = fwl[0].dest;
        auto f1 = fwl[1].dest;
        auto b0 = bwl[0].dest;
        auto b1 = bwl[1].dest;
        if (ws.sg.get_bw_links(f0).size()!=1
            or ws.sg.get_bw_links(f1).size()!=1
               or ws.sg.get_fw_links(b0).size()!=1
                  or ws.sg.get_fw_links(b1).size()!=1)
            continue;
        sgNodeID_t all[4] = {f0, f1, b0, b1};
        bool ok = true;
        for (auto x:all) if (x == n or x == -n or ws.sg.nodes[(x > 0 ? x : -x)].sequence.size() < 199) ok = false;
        for (auto j = 0; j < 3; ++j)
            for (auto i = j + 1; i < 4; ++i)
                if (all[i] == all[j] or all[i] == -all[j])ok = false; //looping node
        if (!ok) continue;

        auto cif0=ws.kci.compute_compression_for_node(f0);
        if (cif0<min_ci or cif0>max_ci) continue;
        auto cif1=ws.kci.compute_compression_for_node(f1);
        if (cif1<min_ci or cif1>max_ci) continue;
        auto cib0=ws.kci.compute_compression_for_node(b0);
        if (cib0<min_ci or cib0>max_ci) continue;
        auto cib1=ws.kci.compute_compression_for_node(b1);
        if (cib1<min_ci or cib1>max_ci) continue;


        auto f0t=ws.linked_read_mappers[0].get_node_tags(f0); if (f0t.size()<3) continue;
        auto f1t=ws.linked_read_mappers[0].get_node_tags(f1); if (f1t.size()<3) continue;
        auto b0t=ws.linked_read_mappers[0].get_node_tags(b0); if (b0t.size()<3) continue;
        auto b1t=ws.linked_read_mappers[0].get_node_tags(b1); if (b1t.size()<3) continue;

        uint64_t aa=intersection_size(b0t,f0t);
        uint64_t bb=intersection_size(b1t,f1t);
        uint64_t ab=intersection_size(b0t,f1t);
        uint64_t ba=intersection_size(b1t,f0t);

        if (aa > 3 and bb > 3 and std::min(aa, bb) > 10 * std::max(ab, ba)) {
            ++aa_count;
            ws.sg.expand_node(n,{{b0},{b1}},{{f0},{f1}});
        } else if (ba > 3 and ab > 3 and
                   std::min(ba, ab) > 10 * std::max(aa, bb)) {
            ++ab_count;
            ws.sg.expand_node(n,{{b0},{b1}},{{f1},{f0}});
        }
        else {
            ++unsolved_count;
        }
    }
    std::cout<<"Repeat expansion summary AA:"<<aa_count<<" AB:"<<ab_count
             <<" Unsolved:"<<unsolved_count<<std::endl;
    return aa_count+ab_count;
}

std::vector<std::pair<sgNodeID_t,sgNodeID_t>> Untangler::get_all_HSPNPs() {
    std::vector<std::pair<sgNodeID_t,sgNodeID_t>> hps;
    std::vector<bool> used(ws.sg.nodes.size(),false);
    //TODO: check the coverages are actually correct?
    const double min_c1=0.75,max_c1=1.25,min_c2=1.5,max_c2=2.5;
    /*
     * the loop always keep the first and the last elements as c=2 collapsed nodes.
     * it starts with a c=2 node, and goes thorugh all bubbles fw, then reverts the subgraph and repeats
     */

#pragma omp parallel for schedule(static, 10)
    for (auto n=1;n<ws.sg.nodes.size();++n){
        std::pair<sgNodeID_t,sgNodeID_t> hap={0,0};
        if (ws.sg.nodes[n].status==sgNodeDeleted) continue;
        auto frontkci=ws.kci.compute_compression_for_node(n);
        if (frontkci>max_c2 or frontkci<min_c2) continue;
        auto m=n;
        //two passes: 0->fw, 1->bw,
        for (auto pass=0; pass<2; ++pass,m=-m) {

            //check bubble going forward ------------
            auto fw_l = ws.sg.get_fw_links(m);
            //fork opening
            if (fw_l.size() != 2) continue;
            hap.first = fw_l[0].dest;
            hap.second = fw_l[1].dest;
            bool used_hspnp;
            #pragma omp critical(find_hspnp_used)
            {
                used_hspnp = (used[(hap.first > 0 ? hap.first : -hap.first)] or
                              used[(hap.second > 0 ? hap.second : -hap.second)]);
                if (!used_hspnp){
                    used[(hap.first>0?hap.first:-hap.first)]=true;
                    used[(hap.second>0?hap.second:-hap.second)]=true;
                }
            }
            if (used_hspnp) continue;
            if (hap.first == n or hap.first == -n or hap.second == n or hap.second == -n or hap.first == hap.second) continue;
            //fork.closing
            auto hap0f = ws.sg.get_fw_links(hap.first);
            auto hap1f = ws.sg.get_fw_links(hap.second);
            if (hap0f.size() != 1 or hap1f.size() != 1 or hap0f[0].dest != hap1f[0].dest) continue;
            auto h0kc = ws.kci.compute_compression_for_node(hap.first);
            if (h0kc > max_c1 or h0kc < min_c1) continue;
            auto h1kc = ws.kci.compute_compression_for_node(hap.second);
            if (h1kc > max_c1 or h1kc < min_c1) continue;
            auto ekc = ws.kci.compute_compression_for_node(hap0f[0].dest);
            if (ekc > max_c2 or ekc < min_c2) continue;

#pragma omp critical(hps)
            {
                hps.push_back(hap);
                if (hps.size()%1000==0) sglib::OutputLog()<<hps.size()<<" haplotype pairs found"<<std::endl;
            }
        }
    }
    sglib::OutputLog()<<hps.size()<<" haplotype pairs found (DONE)"<<std::endl;
    return hps;
}


/**
 * @brief generates tag walk for each HSPNPs, finds common extensions and extends the bubbles
 * @return
 */
uint64_t Untangler::extend_HSPNPs_by_tagwalking() {
    auto const hps=get_all_HSPNPs();
    //std::cout<<"limiting HSPNs to 1000!!! (for test purposes)"<<std::endl;
    //hps.resize(100);
    std::atomic<uint64_t> processing(0);
    std::vector<std::vector<SequenceGraphPath>> new_ppaths;
#pragma omp parallel for
    for (auto i=0;i<hps.size();++i) {
        uint64_t p;
        if ((p=++processing)%100==0) sglib::OutputLog()<<"Procesing HSPNP #"<<p<<std::endl;
        //if (ws.sg.nodes[llabs(hp.first)].sequence.size()<2000 or ws.sg.nodes[llabs(hp.second)].sequence.size()<2000 ) continue;
        TagWalker tw(ws,hps[i]);
        auto ct= tw.remove_crosstalk();
        if (ct>0.01) continue;
        auto wp=tw.walk(.98,.02);
        auto parallelpaths=make_parallel_paths(wp);
        //tw.dump_reads("HPSNP_"+std::to_string(llabs(hp.first))+"_"+std::to_string(llabs(hp.second)));


        //walk_from(hp.first,ws);
        //walk_from(hp.second,ws);
#pragma omp critical (print_paths)
        {

            //std::cout << std::endl << "PATH A (" << wp[0].get_sequence().size() << " bp): ";
            //for (auto n:wp[0].nodes) std::cout << "seq" << llabs(n) << ", ";
            //std::cout << std::endl << "PATH B (" << wp[1].get_sequence().size() << " bp): ";
            //for (auto n:wp[1].nodes) std::cout << "seq" << llabs(n) << ", ";
            //std::cout << std::endl;
            if (parallelpaths[0].nodes.size()<4 or not all_nodes_consumed(parallelpaths)){
                //std::cout << std::endl <<"NO ALL-CONSUMING PARALLEL PATHS"<<std::endl;
            }
            else {
                new_ppaths.emplace_back(parallelpaths);
                //std::cout << "PARALLEL PATHS #"<<new_ppaths.size()-1;
                //std::cout << std::endl << "PARALLEL PATH A (" << parallelpaths[0].get_sequence().size() << " bp): ";
                //for (auto n:parallelpaths[0].nodes) std::cout << "seq" << llabs(n) << ", ";
                //std::cout << std::endl << "PARALLEL PATH B (" << parallelpaths[1].get_sequence().size() << " bp): ";
                //for (auto n:parallelpaths[1].nodes) std::cout << "seq" << llabs(n) << ", ";
                //std::cout << std::endl;
            }
        }

    }
    //get the "best neighbour" for each pp
//    for (auto i1=0;i1<new_ppaths.size();++i1){
//        for (auto i2=i1+1;i2<new_ppaths.size();++i2){
//            auto shared=shared_nodes({new_ppaths[i1],new_ppaths[i2]});
//            if (shared.size()>1) {
//                std::cout<<"Parallel paths #"<<i1<<" and #"<<i2<<" share "<<shared.size()<<" nodes"<<std::endl;
//            }
//        }
//    }
    //execute merging/conflict

    //now join the paths that survived
    std::sort(new_ppaths.begin(),new_ppaths.end(),
              [](const std::vector<SequenceGraphPath> &a,const std::vector<SequenceGraphPath> &b)
              {return a[0].nodes.size()>b[0].nodes.size();});
    std::vector<bool> used(ws.sg.nodes.size());
    for (auto pp:new_ppaths){
        if (pp[0].nodes.size()<4) continue;
        bool skip=false;
        for (auto &p:pp) for (auto n:p.nodes) if (used[llabs(n)]) {skip=true; break;}
        if (skip) continue;
        std::cout << "JOINING PARALLEL PATHS #"<<new_ppaths.size()-1;
        std::cout << std::endl << "PARALLEL PATH A (" << pp[0].get_sequence().size() << " bp): ";
        for (auto n:pp[0].nodes) std::cout << "seq" << llabs(n) << ", ";
        std::cout << std::endl << "PARALLEL PATH B (" << pp[1].get_sequence().size() << " bp): ";
        for (auto n:pp[1].nodes) std::cout << "seq" << llabs(n) << ", ";
        std::cout << std::endl;
        for (auto &p:pp) {
            std::vector<sgNodeID_t> middle_nodes;
            for (auto it=p.nodes.begin()+1;it<p.nodes.end()-1;++it) middle_nodes.push_back(*it);
            for (auto n:middle_nodes) used[llabs(n)]=true;
            ws.sg.join_path(SequenceGraphPath(ws.sg,middle_nodes),true);
        }
    }
    ws.sg.write_to_gfa("graph_after_joining_walks.gfa");
    ws.linked_read_mappers[0].remove_obsolete_mappings();
    ws.linked_read_mappers[0].update_graph_index();
    ws.linked_read_mappers[0].map_reads();
    ws.kci.reindex_graph();
}

/**
 * @brief finds a common source and a common sink, so the paths run parallel through a region
 * @return new paths that start and finish on the common source and sink nodes
 */
std::vector<SequenceGraphPath> Untangler::make_parallel_paths(std::vector<SequenceGraphPath> paths){
    sgNodeID_t source=0,sink=0;
    for (auto na:paths[0].nodes){
        bool missing=false;
        for (auto op=paths.begin()+1;op<paths.end();++op)
            if (std::find(op->nodes.begin(),op->nodes.end(),na)==op->nodes.end()) {
                missing=true;
                break;
            }
        if (missing) continue;
        if (source==0) source=na;
        sink=na;
    }
    std::vector<SequenceGraphPath> pp;
    for (auto p:paths){
        pp.emplace_back(ws.sg);
        bool started=false,finished=false;
        for (auto n:p.nodes){
            if (n==source) started=true;
            if(started and not finished) pp.back().nodes.emplace_back(n);
            if (n==sink) finished=true;
        }
    }
    return pp;
}

/**
 * @brief validates that a parallel path set uses all nodes between the common source and sink
 * @return
 */
bool Untangler::all_nodes_consumed(std::vector<SequenceGraphPath> parallel_paths){
    std::unordered_set<sgNodeID_t> all_nodes;
    for (auto &p:parallel_paths) for (auto &n:p.nodes) {all_nodes.insert(n);all_nodes.insert(-n);}

    for (auto &p:parallel_paths){
        for (auto ix=p.nodes.begin()+1;ix<p.nodes.end()-1;++ix) {
            for (auto fc:ws.sg.get_fw_links(*ix))if (all_nodes.count(fc.dest) == 0) return false;
            for (auto bc:ws.sg.get_bw_links(*ix))if (all_nodes.count(bc.dest) == 0) return false;
        }

    }
    return true;

}

std::vector<sgNodeID_t> Untangler::shared_nodes(std::vector<std::vector<SequenceGraphPath>> parallel_paths){
    std::set<sgNodeID_t> seen_nodes, current_nodes, shared;
    std::vector<sgNodeID_t > r;
    //checks if any of the nodes are present in more than 1 of the PP
    for (auto &pp:parallel_paths){
        current_nodes.clear();
        for (auto &p:pp) for (auto n:p.nodes) current_nodes.insert(llabs(n));
        std::set_intersection(seen_nodes.begin(),seen_nodes.end(),current_nodes.begin(),current_nodes.end(),std::inserter(shared,shared.end()));
        for (auto &n:current_nodes) seen_nodes.insert(n);
    }
    for (auto n:shared) r.emplace_back(n);
    return r;
}

std::vector<SequenceGraphPath> Untangler::combine(std::vector<SequenceGraphPath> parallel_paths1, std::vector<SequenceGraphPath> parallel_paths2){
//    bool fw=false,rev=false;
//    int shared[parallel_paths1.size()][parallel_paths2.size()];
//    int i1=0,i2=0;
//    for (auto & pp1:parallel_paths1) {
//        i2=0;
//        for  (auto & pp2:parallel_paths2){
//            shared[i1][i2]=0;
//            for (auto n1:pp1..nodes) {
//                for (auto n2:pp2.nodes) {
//                    if (n1==n2) {
//                        fw=true;
//                        ++shared[i1][i2];
//                    }
//                    if (n1==-n2) {
//                        rev=true;
//                        ++shared[i1][i2];
//                    }
//                }
//            }
//            ++i2;
//        }
//        ++i1;
//    }
//    if (fw and rev) return {};
//    if (!fw and !rev) return {};
//    if (rev) for (auto &pp:parallel_paths2) pp.reverse(); //now we're forward;
//
//    if (parallel_paths1.size()!=2 or parallel_paths2.size()!=2){
//        std::cout<<"Warning: parallel paths not merging because size>2"<<std::endl;
//        return {};
//    }
//    auto & a1=parallel_paths1[0];
//    auto & b1=parallel_paths1[1];
//    auto & a2=parallel_paths2[0];
//    auto & b2=parallel_paths2[1];
//    if (shared[0][1]>shared[0][0] and shared[1][0]>shared[1][1]){
//        std::swap(a2.nodes,b2.nodes);
//    }
//    //
//    int aoff=-1;
//    for (auto i=0;i<a1.nodes.size();++i) {
//        if (a1.nodes[i] == a2.nodes[0]) {
//            aoff = i;
//            break;
//        }
//    }
//    if (aoff<0) std::swap(a2.nodes,a1.nodes);
//    for (auto i=0;i<a1.nodes.size();++i) {
//        if (a1.nodes[i] == a2.nodes[0]) {
//            aoff = i;
//            break;
//        }
//    }
//    if (aoff<0) return {};
//
//    //pp1.uniqueA intersection pp2.uniqueA? pp2.uniqueB? -pp2.uniqueA? -pp2.uniqueB?
//
//    //set A1 A2, B1, B2 to be the paths that must be combined.
//
//    //for each side
//        //find start of overlap
//        //check overlap ==

}


void Untangler::analise_paths_through_nodes() {
//    std::vector<uint8_t> used(ws.sg.nodes.size(),0);
//    for (auto &p:ws.path_datastores[0].paths){
//        for (auto i=1;i<p.nodes.size()-1;++i) if (used[llabs(p.nodes[i])]<20) ++used[llabs(p.nodes[i])];
//    }
//    uint64_t hall[21];
//    for (auto &x:hall) x=0;
//    uint64_t h399[21];
//    for (auto &x:h399) x=0;
//    uint64_t h402[21];
//    for (auto &x:h402) x=0;
//    uint64_t h1000[21];
//    for (auto &x:h1000) x=0;
//    uint64_t h10000[21];
//    for (auto &x:h10000) x=0;
//    for (auto n=0;n<ws.sg.nodes.size();++n) {
//        ++hall[used[n]];
//        if (ws.sg.nodes[n].sequence.size()>=399) {
//            if (ws.sg.nodes[n].sequence.size()<402) ++h399[used[n]];
//            else if (ws.sg.nodes[n].sequence.size()<1000) ++h402[used[n]];
//            else if (ws.sg.nodes[n].sequence.size()<10000) ++h1000[used[n]];
//            else ++h10000[used[n]];
//        }
//
//    }
//    std::cout<<std::endl<<"=== Paths through -> nodes ==="<<std::endl;
//    for (auto i=0;i<21;++i) std::cout<<i<<", "<<hall[i]<<std::endl;
//    std::cout<<std::endl<<"=== Paths through -> nodes 399bp - 401bp ==="<<std::endl;
//    for (auto i=0;i<21;++i) std::cout<<i<<", "<<h399[i]<<std::endl;
//    std::cout<<std::endl<<"=== Paths through -> nodes 402bp - 999bp ==="<<std::endl;
//    for (auto i=0;i<21;++i) std::cout<<i<<", "<<h402[i]<<std::endl;
//    std::cout<<std::endl<<"=== Paths through -> nodes 1000bp - 9999bp ==="<<std::endl;
//    for (auto i=0;i<21;++i) std::cout<<i<<", "<<h1000[i]<<std::endl;
//    std::cout<<std::endl<<"=== Paths through -> nodes 10000bp+ ==="<<std::endl;
//    for (auto i=0;i<21;++i) std::cout<<i<<", "<<h10000[i]<<std::endl;
//TODO: this is the path-through part
//    struct path_through{
//        sgNodeID_t from,to;
//        uint16_t count;
//        path_through(sgNodeID_t _f,sgNodeID_t _t, uint16_t _c):from(_f),to(_t),count(_c){};
//        const bool operator==(std::pair<sgNodeID_t,sgNodeID_t> a){
//            return (from == a.first and to == a.second) or (from == a.second and to == a.first);
//        }
//    };
//    sglib::OutputLog()<<"Creating path_through structures for each node"<<std::endl;
//    std::vector<std::vector<struct path_through>> paths_through_nodes(ws.sg.nodes.size());
//    std::pair<sgNodeID_t,sgNodeID_t> pfromto;
//    //TODO: add the first and last with origin/end in themselves
//    for (auto &p:ws.path_datastores[0].paths){
//        for (auto i=1;i<p.nodes.size()-1;++i) {
//            auto n=llabs(p.nodes[i]);
//            if (p.nodes[i]>0) { pfromto.first=p.nodes[i-1];pfromto.second=-p.nodes[i+1]; }
//            else { pfromto.first=-p.nodes[i+1]; pfromto.second=p.nodes[i-1];}
//
//            auto pti=std::find(paths_through_nodes[n].begin(),paths_through_nodes[n].end(),pfromto);
//            if (pti==paths_through_nodes[n].end()) paths_through_nodes[n].emplace_back(pfromto.first,pfromto.second,1);
//            else ++(pti->count);
//        }
//    }

    //Canonical repeats.
//    for (sgNodeID_t n=0;n<ws.sg.nodes.size();++n){
//        auto fwl=ws.sg.get_fw_links(n);
//        auto bwl=ws.sg.get_bw_links(n);
//        if (fwl.size()>1 and fwl.size()==bwl.size() ){
//            if (paths_through_nodes[n].size()!=fwl.size()) continue;
//            std::set<sgNodeID_t> in_sol;
//            for (auto &ptn:paths_through_nodes[n]) {
//                in_sol.insert(llabs(ptn.from));
//                in_sol.insert(llabs(ptn.to));
//            }
//            if (in_sol.size()<2*fwl.size()) continue;
//            std::cout<<std::endl<<"Canonical repeat x"<<fwl.size()<<" on node "<<n<<" ("<<ws.sg.nodes[n].sequence.size()<<"bp) solved?"<<std::endl;
//            for (auto &ptn:paths_through_nodes[n]) std::cout<<" "<<ptn.from<<" -> "<<ptn.to<<": "<<ptn.count<<std::endl;
//        }
//
//    }
    std::cout<<"computing KCI for all nodes"<<std::endl;
    ws.kci.compute_all_nodes_kci(1);
    std::ofstream kciof("kci_dump.csv");
    for (auto i=1;i<ws.sg.nodes.size();++i) kciof<<"seq"<<i<<", "<<ws.kci.nodes_depth[i]<<std::endl;
    std::cout<<"DONE!"<<std::endl;
    //CRAP regions
    std::vector<bool> used(ws.sg.nodes.size(),false),aborted(ws.sg.nodes.size(),false);
    std::vector<SequenceSubGraph> craps;
    for (sgNodeID_t n=1;n<ws.sg.nodes.size();++n) {
        if (used[n]) continue;
        if (ws.kci.nodes_depth[n]<1.5 and ws.sg.nodes[n].sequence.size()>500) continue;
        SequenceSubGraph crap(ws.sg);
        std::vector<sgNodeID_t> to_explore={n};
        std::cout<<"Exploring node "<<n<<std::endl;
        while (!to_explore.empty()){
            std::vector<sgNodeID_t> new_to_explore;
            //std::cout<<" Exploring "<<to_explore.size()<< " neighbors"<<std::endl;
            for (auto ne:to_explore){
                //std::cout<<"  Exploring node #"<<ne<<std::endl;
                for (auto neigh:ws.sg.links[ne]){
                    auto x=llabs(neigh.dest);
                    //std::cout<<"  Exploring neighbour #"<<x<<std::endl;
                    if (std::find(crap.nodes.begin(),crap.nodes.end(),x)==crap.nodes.end()) {
                        //std::cout<<"  ADDING to crap!"<<x<<std::endl;
                        crap.nodes.emplace_back(x);
                        if (not (ws.kci.nodes_depth[x]<1.5 and ws.sg.nodes[x].sequence.size()>500)) {
                            //std::cout<<"  ADDING to explore list (size="<<ws.sg.nodes[x].sequence.size()<<")"<<std::endl;
                            new_to_explore.emplace_back(x);
                        }
                    }
                }
            }
            //std::cout<<new_to_explore.size()<< " new neighbors to explore on the next round"<<std::endl;
            to_explore=new_to_explore;
            if (crap.nodes.size()>5000) break;
        }
        if (crap.nodes.size()>10) {
            for (auto &x:crap.nodes) used[x]=true;
            if (crap.nodes.size()<=5000) {
                craps.emplace_back(crap);
                std::cout << std::endl << "==== NEW CRAP ====" << std::endl;
                for (auto x:crap.nodes) std::cout << "seq" << x << ", ";
                std::cout << std::endl;
            }
            if (crap.nodes.size()>5000) {
                for (auto x:crap.nodes) aborted[x]=true;
            }
        }
    }
    uint64_t abt=0;
    for (auto a:aborted) if (a) ++abt;
    std::cout<<"There were "<<abt<<" nodes in aborted crap components"<<std::endl;


}

std::vector<std::pair<sgNodeID_t,sgNodeID_t>> Untangler::find_bubbles(uint32_t min_size,uint32_t max_size) {
    std::vector<std::pair<sgNodeID_t,sgNodeID_t>> r;
    std::vector<bool> used(ws.sg.nodes.size(),false);
    sgNodeID_t n1,n2;
    size_t s1,s2;
    for (n1=1;n1<ws.sg.nodes.size();++n1){
        if (used[n1]) continue;
        //get "topologically correct" bubble: prev -> [n1 | n2] -> next
        s1=ws.sg.nodes[n1].sequence.size();
        if (s1<min_size or s1>max_size) continue;
        auto fwl=ws.sg.get_fw_links(n1);
        if (fwl.size()!=1) continue;
        auto next=fwl[0].dest;
        auto bwl=ws.sg.get_bw_links(n1);
        if (bwl.size()!=1) continue;
        auto prev=bwl[0].dest;
        auto parln=ws.sg.get_bw_links(next);
        if (parln.size()!=2) continue;
        auto parlp=ws.sg.get_fw_links(-prev);
        if (parlp.size()!=2) continue;
        if (parlp[0].dest!=n1) n2=parlp[0].dest;
        else n2=parlp[1].dest;
        if (n2!=-parln[0].dest and n2!=-parln[1].dest) continue;
        s2=ws.sg.nodes[n2].sequence.size();
        if (s2<min_size or s2>max_size) continue;
        used[n1]=true;
        used[n2]=true;
        r.emplace_back(n1,n2);
    }
    return r;
}

std::vector<std::pair<sgNodeID_t,sgNodeID_t>> Untangler::solve_bubbly_paths() {
    //find bubbly paths
    auto bps=ws.sg.get_all_bubbly_subgraphs();
    std::vector<SequenceSubGraph> kobps,sbps;
    uint64_t allnodes=0,bnodes=0,allnodes_size=0,bnodes_size=0;
    for (auto &bp:bps){
        /*if (bp.nodes.size()>40) {
            std::cout<<"LARGE bubbly path: ";
            for (auto n:bp.nodes) std::cout<<"seq"<<labs(n)<<",";
            std::cout<<std::endl;

        }*/
        //check every third node for kci
        allnodes+=bp.nodes.size();
        allnodes_size+=bp.total_size();
        int kci_ok=0,kci_fail=0;
        for (auto i=0;i<bp.nodes.size();i+=3) {
            auto dup_kci=ws.kci.compute_compression_for_node(bp.nodes[i],1);
            if (dup_kci<1.5 or dup_kci>2.5) ++kci_fail;
            else ++kci_ok;
        }
        if (kci_fail*3<kci_ok) {
            kobps.push_back(bp);
            bnodes+=bp.nodes.size();
            bnodes_size+=bp.total_size();
        }
    }
    std::cout<<kobps.size()<<" / "<< bps.size() << " bubly subgraphs ( "<<bnodes<<" / "<<allnodes<<" nodes) ( "<<bnodes_size<<" / "<<allnodes_size<<" bp) pass kci condition on duplicated nodes"<<std::endl;
    uint64_t solved=0, solved_nodes=0, solved_size=0, untagged=0,untagged_nodes=0, untagged_size=0 ,ambiguous=0,ambiguous_nodes=0, ambiguous_size=0;
    std::vector<std::pair<std::vector<sgNodeID_t>,std::vector<sgNodeID_t>>> solved_haps;
    int min_tags=10;
    for (auto &bp:kobps){
        //init the haplotypes anchoring all nodes to previous tags
        //std::cout<<"A";
        std::vector<sgNodeID_t> hap1, hap2;
        std::set<bsg10xTag> tags1,tags2;
        hap1.push_back(bp.nodes[1]);
        tags1=ws.linked_read_mappers[0].get_node_tags(bp.nodes[1]);
        hap2.push_back(bp.nodes[2]);
        tags2=ws.linked_read_mappers[0].get_node_tags(bp.nodes[2]);
        if (tags1.size()<min_tags or tags2.size()<min_tags) {++untagged; untagged_nodes+=bp.nodes.size(); untagged_size+=bp.total_size(); continue;};
        bool utb=false;
        for (auto i=3;i<bp.nodes.size()-1;i+=3){
            std::set<bsg10xTag> both;
            std::set_intersection(tags1.begin(),tags1.end(),tags2.begin(),tags2.end(),std::inserter(both,both.end()));
            for (auto b:both){
                tags1.erase(b);
                tags2.erase(b);
            }
            auto A=bp.nodes[i+1];
            auto B=bp.nodes[i+2];
            auto tagsA=ws.linked_read_mappers[0].get_node_tags(A);
            auto tagsB=ws.linked_read_mappers[0].get_node_tags(B);
            if (tagsA.size()<min_tags or tagsB.size()<min_tags) {++untagged;untagged_nodes+=bp.nodes.size(); untagged_size+=bp.total_size(); utb=true; break;};
            auto a1=intersection_size(tagsA,tags1);
            auto a2=intersection_size(tagsA,tags2);
            auto b1=intersection_size(tagsB,tags1);
            auto b2=intersection_size(tagsB,tags2);
            if (a1>=min_tags and b1<a1*2/min_tags and b2>=min_tags and a2<b2*2/min_tags){
                //std::cout<<"A";
                hap1.push_back(A);
                tags1.insert(tagsA.begin(),tagsA.end());
                hap2.push_back(B);
                tags2.insert(tagsB.begin(),tagsB.end());
            }
            else if (b1>=min_tags and a1<b1*2/min_tags and a2>=min_tags and b2<a2*2/min_tags){
                //std::cout<<"B";
                hap2.push_back(A);
                tags2.insert(tagsA.begin(),tagsA.end());
                hap1.push_back(B);
                tags1.insert(tagsB.begin(),tagsB.end());
            } else {
                //std::cout<<" ["<<a1<<","<<a2<<","<<b1<<","<<b2<<"] FAIL!";
                break;
            }
        }
        //std::cout<<std::endl;
        if (utb) continue;
        if (hap1.size()==bp.nodes.size()/3) {
            ++solved;
            sbps.push_back(bp);
            solved_haps.push_back({hap1,hap2});
            solved_nodes+=bp.nodes.size();
            solved_size+=bp.total_size();
        }
        else{
            ++ambiguous;
            ambiguous_nodes+=bp.nodes.size();
            ambiguous_size+=bp.total_size();
        }

        //
        //for (auto)

    }
    std::cout<<untagged << " / " << kobps.size()<<" bubly subgraphs had <"<<min_tags<<" tags in a unique node (total "<<untagged_nodes<<" nodes / "<<untagged_size<<" bp)"<<std::endl;
    std::cout<<ambiguous << " / " << kobps.size()<<" bubly subgraphs ambiguous (total "<<ambiguous_nodes<<" nodes) / "<<ambiguous_size<<" bp)"<<std::endl;
    std::cout<<solved << " / " << kobps.size()<<" bubly subgraphs solved (total "<<solved_nodes<<" nodes) / "<<solved_size<<" bp)"<<std::endl;

    solved=0;
    //expand each repeat
    for (auto i=0;i<sbps.size();++i){
        //std::cout<<"Solving region: ";
        //for (auto n:sbps[i].nodes) std::cout<<n<<" ";
        //std::cout<<std::endl;
        //std::cout<<"Hap1: ";
        //for (auto n:solved_haps[i].first) std::cout<<n<<" ";
        //std::cout<<std::endl;
        //std::cout<<"Hap2: ";
        //for (auto n:solved_haps[i].second) std::cout<<n<<" ";
        //std::cout<<std::endl;
        for (auto l=0;l<solved_haps[i].first.size()-1;++l){
            auto r=sbps[i].nodes[3*(l+1)];
            auto a1=solved_haps[i].first[l];
            auto a2=solved_haps[i].first[l+1];
            auto b1=solved_haps[i].second[l];
            auto b2=solved_haps[i].second[l+1];
            //std::cout<<"Expanding "<<a1<<" -> "<<r<<" -> "<<a2<<"   |    "<<b1<<" -> "<<r<<" -> "<<b2<<std::endl;
            ws.sg.expand_node(r,{{a1},{b1}},{{a2},{b2}});
            ++solved;
        }

    }
    std::cout<<solved << " duplications expanded"<<std::endl;
    //done!
    return {};
}

void Untangler::pop_errors_by_ci_and_paths() {
    sglib::OutputLog()<<"Popping errors..."<<std::endl;
    auto bubbles=find_bubbles(200, 450);
    sglib::OutputLog()<<"Analysing "<<bubbles.size()<<" small bubbles for coverage"<<std::endl;
    for (auto bp:bubbles){
        auto ci1=ws.kci.compute_compression_for_node(bp.first);
        auto ci2=ws.kci.compute_compression_for_node(bp.second);
        if (ci1>.7 and ci2<.1) std::cout<<"node "<<bp.second<<" has only "<<ci2<<" coverage and "<<bp.first<<" has "<<ci1<<std::endl;
        if (ci2>.7 and ci1<.1) std::cout<<"node "<<bp.first<<" has only "<<ci1<<" coverage and "<<bp.second<<" has "<<ci2<<std::endl;

    }
}

/**
 * @brief grabs all "long" haplotype-specific nodes, uses tags to find neighbours.
 * @param min_size
 * @param min_ci
 * @param max_ci
 * @return
 */
std::vector<std::vector<std::pair<sgNodeID_t,uint32_t>>> Untangler::find_tag_neighbours(uint32_t min_size, float min_ci, float max_ci) {


    sglib::OutputLog()<<"Selecting nodes..."<<std::endl;
    std::vector<std::vector<std::pair<sgNodeID_t,uint32_t>>> neighbours;
    std::vector<std::set<bsg10xTag>> node_tags;
    neighbours.resize(ws.sg.nodes.size());
    node_tags.resize(ws.sg.nodes.size());
    uint64_t total_bp=0;
    auto nodes=ws.select_from_all_nodes(min_size,1000000,20,200000, min_ci, max_ci);
    sglib::OutputLog()<<"Populating node tags..."<<std::endl;
    for (auto n:nodes) {
        total_bp+=ws.sg.nodes[n].sequence.size();
        for (auto t:ws.linked_read_mappers[0].get_node_tags(n)) node_tags[n].insert(t);
    }
    sglib::OutputLog()<<nodes.size()<<" selected totalling "<<total_bp<<"bp "<<std::endl;
    sglib::OutputLog()<<"Computing shared tags"<<std::endl;
#pragma omp parallel for shared(neighbours) schedule(static,50)
    for (auto i1=0; i1<nodes.size(); ++i1){
        auto n1=nodes[i1];
        for (auto i2=i1+1;i2<nodes.size();i2++){
            auto n2=nodes[i2];
            uint32_t shared=intersection_size(node_tags[n1],node_tags[n2]);
            if (shared>=4) {
#pragma omp critical
                {
                    neighbours[n1].emplace_back(n2,shared);
                    neighbours[n2].emplace_back(n1,shared);
                }
            }
        }
    }
    sglib::OutputLog()<<"Sorting shared tags"<<std::endl;
    uint64_t with_neighbours=0;
    for (auto &nn:neighbours){
        if (not nn.empty()) ++with_neighbours;
        std::sort(nn.begin(),nn.end(),[]( const std::pair<sgNodeID_t,uint32_t> &a,
                                          const std::pair<sgNodeID_t,uint32_t> &b) { return a.second>b.second; });
    }
    sglib::OutputLog()<<with_neighbours<<" nodes with neighbours"<<std::endl;
    return neighbours;
}

void Untangler::connect_neighbours() {
    //first find all nodes' neighbours
    auto tagneighbours=find_tag_neighbours(5000,.75,1.25);
    TagWalker tw(ws,{});
    for (auto n=1;n<ws.sg.nodes.size();++n) {
        //explore from a node til hitting a neighbour, check if another selected node is on the way
        std::set<sgNodeID_t> ntn;
        for (auto nd:tagneighbours[n]) ntn.emplace(nd.first);
        if (ntn.empty()) continue;
        //try to skate from a neighbour to another
        std::cout<<std::endl<<"Finding forward neighbours for "<<n<<std::endl;
        auto ndfs=ws.sg.get_distances_to(n,ntn,100000);
        for (auto nd:ndfs){
            std::cout<<"distance to "<<nd.first<<": "<<nd.second<<std::endl;
        }
        std::cout<<std::endl<<"Finding backward neighbours for "<<n<<std::endl;
        auto ndbs=ws.sg.get_distances_to(-n,ntn,100000);
        for (auto nd:ndbs){
            std::cout<<"distance to "<<nd.first<<": "<<nd.second<<std::endl;
        }
        /*//if successfull, create a copy of the skated path and disconnect/connect appropriately
        //TODO: check if a node with different neighbours is here
        auto walkf=tw.walk_between(n,nd[0].first);
        if (not walkf.nodes.empty()) {
            sg.create_path_through(walkf);
        }
        auto nd=ws.sg.get_distances_to(-n,ntn,100000);
        //if successfull, create a copy of the skated path and disconnect/connect appropriately
        //TODO: check if a node with different neighbours is here
        auto walkb=tw.walk_between(n,nd[0].first);
        if (not walkb.nodes.empty()) {
            sg.create_path_through(walkb);
        }*/

    }

}