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
        std::set<bsg10xTag> b0tags, b1tags, f0tags, f1tags;

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

        std::set<bsg10xTag> shared1,shared2;

        for (auto t:f0tags) if (f1tags.count(t)>0) {shared1.insert(t);};
        for (auto t:shared1){f0tags.erase(t);f1tags.erase(t);};
        for (auto t:b0tags) if (b1tags.count(t)>0) {shared2.insert(t);};
        for (auto t:shared2) {b0tags.erase(t);b1tags.erase(t);};

        std::set<bsg10xTag> aa, bb, ba, ab;
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

uint64_t Untangler::expand_canonical_repeats_by_tags(float min_ci, float max_ci, int min_tags) {
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

        auto f0t=ws.linked_read_mappers[0].get_node_tags(f0); if (f0t.size()<min_tags) continue;
        auto f1t=ws.linked_read_mappers[0].get_node_tags(f1); if (f1t.size()<min_tags) continue;
        auto b0t=ws.linked_read_mappers[0].get_node_tags(b0); if (b0t.size()<min_tags) continue;
        auto b1t=ws.linked_read_mappers[0].get_node_tags(b1); if (b1t.size()<min_tags) continue;

        auto cif0=ws.kci.compute_compression_for_node(f0,1);
        if (cif0<min_ci or cif0>max_ci) continue;
        auto cif1=ws.kci.compute_compression_for_node(f1,1);
        if (cif1<min_ci or cif1>max_ci) continue;
        auto cib0=ws.kci.compute_compression_for_node(b0,1);
        if (cib0<min_ci or cib0>max_ci) continue;
        auto cib1=ws.kci.compute_compression_for_node(b1,1);
        if (cib1<min_ci or cib1>max_ci) continue;

        uint64_t aa=intersection_size(b0t,f0t);
        uint64_t bb=intersection_size(b1t,f1t);
        uint64_t ab=intersection_size(b0t,f1t);
        uint64_t ba=intersection_size(b1t,f0t);

        if (aa > min_tags and bb > min_tags and std::min(aa, bb) > 10 * std::max(ab, ba)) {
            ++aa_count;
            ws.sg.expand_node(n,{{b0},{b1}},{{f0},{f1}});
        } else if (ba > min_tags and ab > min_tags and std::min(ba, ab) > 10 * std::max(aa, bb)) {
            ++ab_count;
            ws.sg.expand_node(n,{{b0},{b1}},{{f1},{f0}});
        }
        else {
            ++unsolved_count;
        }
    }
    std::cout<<"Repeat expansion summary AA:"<<aa_count<<" AB:"<<ab_count <<" Unsolved:"<<unsolved_count<<std::endl;
    return aa_count+ab_count;
}

std::vector<std::pair<sgNodeID_t,sgNodeID_t>> Untangler::get_all_HSPNPs() {
    std::vector<std::pair<sgNodeID_t,sgNodeID_t>> hps;
    std::vector<bool> used(ws.sg.nodes.size(),false);
    //TODO: check the coverages are actually correct?
    const double min_c1=0.5,max_c1=1.5,min_c2=1.25,max_c2=3;
    /*
     * the loop always keep the first and the last elements as c=2 collapsed nodes.
     * it starts with a c=2 node, and goes thorugh all bubbles fw, then reverts the subgraph and repeats
     */

#pragma omp parallel for schedule(static, 10)
    for (auto n=1;n<ws.sg.nodes.size();++n){
        std::pair<sgNodeID_t,sgNodeID_t> hap={0,0};
        if (ws.sg.nodes[n].status==sgNodeDeleted) continue;
        //auto frontkci=ws.kci.compute_compression_for_node(n);
        //if (frontkci>max_c2 or frontkci<min_c2) continue;
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
            //auto ekc = ws.kci.compute_compression_for_node(hap0f[0].dest);
            //if (ekc > max_c2 or ekc < min_c2) continue;

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
    ws.sg.create_index();
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




void Untangler::analise_paths_through_nodes() {
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

/**
 * @brief solves a single bubbly path, returns the paths for the 2 new nodes or empty paths if unsolved
 * @param bp
 * @return
 */
std::pair<SequenceGraphPath,SequenceGraphPath> Untangler::solve_bubbly_path(const SequenceSubGraph &bp, bool &no_tags) {
    int min_tags=20; no_tags=false;
    std::vector<std::set<bsg10xTag>> node_tags;

    for (auto i=0;i<bp.nodes.size();++i) {
        auto n=bp.nodes[i];
        node_tags.emplace_back(ws.linked_read_mappers[0].get_node_tags(n));
        if (i%3>0 and node_tags.back().size()<min_tags) {
            no_tags=true;
            return {SequenceGraphPath(ws.sg),SequenceGraphPath(ws.sg)};
        };//check failed: min-tags
    }
    //init the haplotypes anchoring all nodes to previous tags
    //std::cout<<"A";
    std::vector<sgNodeID_t> hap1, hap2;
    std::set<bsg10xTag> tags1,tags2;
    hap1.push_back(bp.nodes[1]);
    tags1=node_tags[1];
    hap2.push_back(bp.nodes[2]);
    tags2=node_tags[2];
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
        auto a1=intersection_size(tagsA,tags1);
        auto a2=intersection_size(tagsA,tags2);
        auto b1=intersection_size(tagsB,tags1);
        auto b2=intersection_size(tagsB,tags2);
        if (a1>=min_tags and b1<a1*2/min_tags and b2>=min_tags and a2<b2*2/min_tags){
            hap1.push_back(A);
            tags1.insert(tagsA.begin(),tagsA.end());
            hap2.push_back(B);
            tags2.insert(tagsB.begin(),tagsB.end());
        }
        else if (b1>=min_tags and a1<b1*2/min_tags and a2>=min_tags and b2<a2*2/min_tags){
            hap2.push_back(A);
            tags2.insert(tagsA.begin(),tagsA.end());
            hap1.push_back(B);
            tags1.insert(tagsB.begin(),tagsB.end());
        } else {
            break;
        }
    }
    //std::cout<<std::endl;
    if (hap1.size()==bp.nodes.size()/3) {
        SequenceGraphPath p1(ws.sg),p2(ws.sg);
        for (auto i=0;i<hap1.size();++i){
            if (i>0) p1.nodes.push_back(bp.nodes[3*i]);
            p1.nodes.push_back(hap1[i]);
            if (i>0) p2.nodes.push_back(bp.nodes[3*i]);
            p2.nodes.push_back(hap2[i]);
        }
        return {p1,p2};
    }
    return {SequenceGraphPath(ws.sg),SequenceGraphPath(ws.sg)};
}

/**
 * @brief solves a single bubbly path, returns the paths for the 2 new nodes or empty paths if unsolved
 * can return many pairs of paths if solving partial parts.
 * @param bp
 * @return
 */
std::vector<std::pair<SequenceGraphPath,SequenceGraphPath>> Untangler::solve_bubbly_path_2(const SequenceSubGraph &bp) {

    std::cout<<std::endl<<"solve_bubbly_path_2 started for "<<bp.nodes.size()<<" nodes ("<<bp.nodes.size()/3<<" bubbles)"
             <<"between "<<bp.nodes.front()<<" and "<<bp.nodes.back()<<std::endl;
    std::cout<<"Nodes in bubbly path:";
    for(auto n:bp.nodes) std::cout<<" "<<n;
    std::cout<<std::endl;

    //CONSTANTS
    int specific_to_shared_tags=100;
    size_t min_score=10;
    size_t max_penalty=3;

    std::set<std::pair<sgNodeID_t, sgNodeID_t>> solution;
    std::map<sgNodeID_t,std::set<bsg10xTag>> node_tags;
    std::set<std::pair<sgNodeID_t, sgNodeID_t>> pnps; //Parallel Node Pairs
    //non-positional best beighbours first: create a set of all PNPs, pick the starting one (highest with up to .5%shared tags
    //TODO:prioritise clean start (shared=0)
    size_t max_min_tags=0;
    std::pair<sgNodeID_t, sgNodeID_t> best_start;
    for (auto i=0;i<bp.nodes.size()-1;i+=3){
        auto A=bp.nodes[i+1];
        auto B=bp.nodes[i+2];
        node_tags[A]=ws.linked_read_mappers[0].get_node_tags(A);
        node_tags[B]=ws.linked_read_mappers[0].get_node_tags(B);
        pnps.insert({A,B});
        auto shared_count=intersection_size(node_tags[A],node_tags[B]);
        if (shared_count*specific_to_shared_tags<node_tags[A].size() and shared_count*specific_to_shared_tags<node_tags[B].size() and max_min_tags<std::min(node_tags[A].size(),node_tags[B].size())){
            best_start.first=A;best_start.second=B;
            max_min_tags=std::min(node_tags[A].size(),node_tags[B].size());
        }
        std::cout<<" bubble tag analysis: A="<<A<<": "<<node_tags[A].size()<<" B="<<B<<": "<<node_tags[B].size()<<"  shared: "<<shared_count<<std::endl;
    }


    if (max_min_tags==0){
        std::cout<<"can't even find where to start, so giving up"<<std::endl;
        return {};
    }

    std::cout<<"starting from PNP "<<best_start.first<<" "<<best_start.second<<std::endl;
    auto tags1=node_tags[best_start.first];
    auto tags2=node_tags[best_start.second];
    auto total_penalty=intersection_size(tags1,tags2);
    auto todo_pnps=pnps;
    todo_pnps.erase(best_start);
    solution.insert(best_start);
    // while possible grab the next lowest-conflict/highest score PNP
    while (not todo_pnps.empty()){
        size_t min_pen=max_penalty;
        size_t max_score=0;
        std::pair<sgNodeID_t, sgNodeID_t> best_next;
        bool flip=false;
        for (auto p:todo_pnps){
            auto a1=intersection_size(node_tags[p.first],tags1);
            auto a2=intersection_size(node_tags[p.first],tags2);
            auto b1=intersection_size(node_tags[p.second],tags1);
            auto b2=intersection_size(node_tags[p.second],tags2);
            //std::cout<<"evaluating PNP ("<<p.first<<", "<<p.second<<") : a1="<<a1<<" a2="<<a2<<" b1="<<b1<<" b2="<<std::endl;
            //TODO: discount already-shared tags from penalty;
            //dir?
            if (a1>min_score and b2> min_score and (a2+b1<min_pen or ( a2+b1==min_pen and a1+b2>max_score))) {
                best_next.first=p.first;best_next.second=p.second;flip=false;
                min_pen=a2+b1;
                max_score=a1+b2;
            }
            //flip?
            else if (a2>min_score and b1> min_score and (a1+b2<min_pen or( a1+b2==min_pen and a2+b1>max_score))) {
                best_next.first=p.first;best_next.second=p.second;flip=true;
                min_pen=a1+b2;
                max_score=a2+b1;
            }
        }
        if (max_score==0) break;
        todo_pnps.erase(best_next);
        if (flip) { auto s=best_next.second;best_next.second=best_next.first;best_next.first=s; }
        solution.insert(best_next);
        for (auto t:node_tags[best_next.first]) tags1.insert(t);
        for (auto t:node_tags[best_next.second]) tags2.insert(t);
    }
    // rescue PNPs without tags for a node (contemplate both error AA BB and non-mappable node walk)
    std::unordered_set<uint64_t> tagkmers1,tagkmers2;

    if (false/*not todo_pnps.empty()*/) {
        std::cout<<"Solution is incomplete, attempting PNP rescue for untagged nodes and sequencing errors"<<std::endl;
        std::set<std::pair<sgNodeID_t ,sgNodeID_t >> rescued;
        StringKMerFactory kfa(31);
        StringKMerFactory kfb(31);
        for (auto p:todo_pnps) {
            auto a1=intersection_size(node_tags[p.first],tags1);
            auto a2=intersection_size(node_tags[p.first],tags2);
            auto b1=intersection_size(node_tags[p.second],tags1);
            auto b2=intersection_size(node_tags[p.second],tags2);
            std::cout<<"Checking PNP ( "<<p.first<<", "<<p.second<<" ): a1="<<a1<<" a2="<<a2<<" b1="<<b1<<" b2="<<b2<<std::endl;
            if (a1>10 and a2>10 and b1<3 and b2<3) {
                std::cout<<"Node "<<p.second<<" seems to be an error!"<<std::endl;
                solution.insert({p.first,p.first});
                rescued.insert(p);
                continue;
            }
            if (b1>10 and b2>10 and a1<3 and a2<3) {
                std::cout<<"Node "<<p.first<<" seems to be an error!"<<std::endl;
                solution.insert({p.second,p.second});
                rescued.insert(p);
                continue;
            }
            //now check for untagged (i.e. kmer coverage)
            if (tagkmers1.empty()){
                std::cout<<"Generating tag kmers for current solution as this is the first walk-in situation"<<std::endl;
                BufferedTagKmerizer btk(ws.linked_read_datastores[0],31,200000,1000);
                std::set<bsg10xTag> exctags1,exctags2;
                for (auto t:tags1) if (tags2.count(t)==0) exctags1.insert(t);
                for (auto t:tags2) if (tags1.count(t)==0) exctags2.insert(t);
                tagkmers1=btk.get_tags_kmers(6, exctags1);
                tagkmers2=btk.get_tags_kmers(6, exctags2);
            }
            std::vector<uint64_t> ka,kb;
            kfa.create_kmers(ws.sg.nodes[llabs(p.first)].sequence, ka);
            kfb.create_kmers(ws.sg.nodes[llabs(p.second)].sequence, kb);
            uint64_t uncovered_a1=0,uncovered_a2=0,uncovered_b1=0,uncovered_b2=0;
            for (auto x:ka) {
                if (tagkmers1.count(x)==0) ++uncovered_a1;
                if (tagkmers2.count(x)==0) ++uncovered_a2;
            }
            for (auto x:kb) {
                if (tagkmers1.count(x)==0) ++uncovered_b1;
                if (tagkmers2.count(x)==0) ++uncovered_b2;
            }
            std::cout<<"kmers not covered by haplotype-specific tags: a~1="<<uncovered_a1<<" a~2="<<uncovered_a2<<" b~1="<<uncovered_b1<<" b~2="<<uncovered_b2<<std::endl;
            if (uncovered_a1==0 and uncovered_b1>0 and uncovered_b2==0 and uncovered_a2>0) {
                std::cout<<"PNP ( "<<p.first<<", "<<p.second<<" ) walked-through as AB!"<<std::endl;
                solution.insert({p.first,p.second});
                rescued.insert(p);
                continue;
            }
            if (uncovered_a2==0 and uncovered_b2>0 and uncovered_b1==0 and uncovered_a1>0) {
                std::cout<<"PNP ( "<<p.first<<", "<<p.second<<" ) walked-through as BA!"<<std::endl;
                solution.insert({p.second,p.first});
                rescued.insert(p);
                continue;
            }
        }
        for (auto p:rescued) todo_pnps.erase(p);
    }
    std::cout<<"Hap1: ";
    for (auto p:solution) std::cout<<"seq"<<llabs(p.first)<<", ";
    std::cout<<std::endl<<"Hap2: ";
    for (auto p:solution) std::cout<<"seq"<<llabs(p.second)<<", ";
    std::cout<<std::endl;
    std::vector<std::pair<SequenceGraphPath,SequenceGraphPath>> paths;
    if (todo_pnps.empty()) {
        std::cout << "FULL SOLUTION found!" << std::endl;
    } else {
        std::cout<<"SOLUTION only contains "<<solution.size()<<"/"<<pnps.size()<<" PNPs"<<std::endl;
        //TODO support partial solutions! (i.e. re-run with broken bubbly path?)
        std::cout<<"Generating partial solutions"<<std::endl;
    }
    std::set<sgNodeID_t> h1nodes,h2nodes;
    for (auto p:solution) {
        h1nodes.insert(p.first);
        h2nodes.insert(p.second);
    }
    std::vector<sgNodeID_t> currh1,currh2;

    for (auto i=0;i<bp.nodes.size();++i){
        if (i%3==0 ){
            if (not currh1.empty()) {
                if (i < bp.nodes.size() - 3 and
                    (h1nodes.count(bp.nodes[i + 1]) > 0 or h2nodes.count(bp.nodes[i + 1]) > 0) and
                    (h1nodes.count(bp.nodes[i + 2]) > 0 or h2nodes.count(bp.nodes[i + 2]) > 0)) {
                    currh1.emplace_back(bp.nodes[i]);
                    currh2.emplace_back(bp.nodes[i]);
                } else {
                    if (currh1.size() > 2) {
                        std::cout << "Solution #" << paths.size() + 1 << ".a:";
                        for (auto n:currh1) std::cout << " " << n;
                        std::cout << std::endl << "Solution #" << paths.size() + 1 << ".b:";
                        for (auto n:currh2) std::cout << " " << n;
                        std::cout << std::endl;
                        paths.push_back({SequenceGraphPath(ws.sg, currh1), SequenceGraphPath(ws.sg, currh2)});
                    }
                    currh1.clear();
                    currh2.clear();
                }
            }
        } else {
            if (h1nodes.count(bp.nodes[i])>0) currh1.emplace_back(bp.nodes[i]);
            if (h2nodes.count(bp.nodes[i])>0) currh2.emplace_back(bp.nodes[i]);
        }
    }
    //remove last node (heterozygous, we don't want to stupidly extend haplotypes



    return paths;
}

std::vector<std::pair<sgNodeID_t,sgNodeID_t>> Untangler::solve_bubbly_paths() {
    //TODO: 3-part structure:
    // 1) find, report validate kci, report
    // 2) For each: solve
    // 3) report

    //find bubbly paths
    auto bps=ws.sg.get_all_bubbly_subgraphs();
    sglib::OutputLog()<<"--- INITIAL bubbly paths ---"<<std::endl;
    ws.sg.print_bubbly_subgraph_stats(bps);
    //TODO: write a more sophisticated kci check
    std::vector<SequenceSubGraph> kobps;
    for (auto &bp:bps){
        int kci_ok=0,kci_fail=0;
        for (auto i=0;i<bp.nodes.size();i+=3) {
            auto dup_kci=ws.kci.compute_compression_for_node(bp.nodes[i],1);
            if (dup_kci<1.25 or dup_kci>2.5) ++kci_fail;
            else ++kci_ok;
        }
        if (kci_fail*3<kci_ok) {
            kobps.push_back(bp);
        }
    }
    sglib::OutputLog()<<"--- KCI OK bubbly paths ---"<<std::endl;
    ws.sg.print_bubbly_subgraph_stats(kobps);



    //uint64_t solved=0, solved_nodes=0, solved_size=0, untagged=0,untagged_nodes=0, untagged_size=0 ,ambiguous=0,ambiguous_nodes=0, ambiguous_size=0;
    //std::vector<std::pair<std::vector<sgNodeID_t>,std::vector<sgNodeID_t>>> solved_haps;
    //int min_tags=20;
    /*{
        uint64_t solved = 0, unsolved = 0, untagged = 0;
        std::vector<SequenceSubGraph> solbubs;
        for (auto &bp:kobps) {
            bool no_tags;
            auto sol = solve_bubbly_path(bp, no_tags);
            if (no_tags) ++untagged;
            if (sol.first.nodes.size() == 0) ++unsolved;
            else {
                ++solved;
                solbubs.push_back(bp);
            }

        }
        sglib::OutputLog() << "OLD solver: "<< solved << " bubbly paths solved, " << unsolved << " unsolved, with " << untagged << " lacking tags" << std::endl;
        ws.sg.print_bubbly_subgraph_stats(solbubs);
        //done!
    }*/
    {
        uint64_t solved = 0, unsolved = 0;
        std::vector<SequenceSubGraph> solbubs;
        std::vector<SequenceGraphPath> paths;
        for (auto &bp:kobps) {
            auto sol = solve_bubbly_path_2(bp);

            if (sol.size() == 0) ++unsolved;
            else {
                ++solved;
                solbubs.push_back(bp);
                for (auto pp:sol) {
                    paths.push_back(pp.first);
                    paths.push_back(pp.second);
                }

            }

        }
        sglib::OutputLog() << "New solver: "<< solved << " bubbly paths solved, " << unsolved << " unsolved" << std::endl;
        ws.sg.print_bubbly_subgraph_stats(solbubs);
        std::set<sgNodeID_t> used_nodes;
        for (auto p:paths) {
            ws.sg.join_path(p,false);
            for (auto n:p.nodes) used_nodes.insert(n);
        }
        for (auto n:used_nodes) ws.sg.remove_node(n);
        //done!
    }
    return {};
}

void Untangler::pop_errors_by_ci_and_paths(uint32_t min_size, uint32_t max_size) {
    sglib::OutputLog()<<"Popping errors..."<<std::endl;
    auto bubbles=find_bubbles(min_size, max_size);
    sglib::OutputLog()<<"Analysing "<<bubbles.size()<<" small bubbles for coverage"<<std::endl;
    std::vector<sgNodeID_t> to_delete;
    //std::ofstream bubblesf("bubbles_detail.csv");
    //bubblesf<<"prev,b1,b2,next,ci_prev,ci1,ci2,ci_next"<<std::endl;
    for (auto bp:bubbles){
        auto ci1=ws.kci.compute_compression_for_node(bp.first,1);
        auto ci2=ws.kci.compute_compression_for_node(bp.second,1);
        auto prev=ws.sg.get_bw_links(bp.first)[0].dest;
        auto next=ws.sg.get_fw_links(bp.first)[0].dest;
        auto cip=ws.kci.compute_compression_for_node(prev);
        auto cin=ws.kci.compute_compression_for_node(next);
        //bubblesf<<prev<<", "<<bp.first<<", "<<bp.second<<", "<<next<<", "<<cip<<", "<<ci1<<", "<<ci2<<", "<<cin<<std::endl;
        if (cip<1.3 and cin<1.3) {
            if (ci1 > .8 and ci2 < .3) {
                //std::cout << "node " << bp.second << " has only " << ci2 << " coverage and " << bp.first << " has "
                //          << ci1 << std::endl;
                to_delete.push_back(llabs(bp.second));
            }
            if (ci2 > .8 and ci1 < .3) {
                //std::cout << "node " << bp.first << " has only " << ci1 << " coverage and " << bp.second << " has "
                //          << ci2 << std::endl;
                to_delete.push_back(llabs(bp.first));
            }
        }
    }
    //std::cout<<"Deleting "<<to_delete.size()<<" nodes as errors"<<std::endl;
    for (auto &pb:to_delete) ws.sg.remove_node(pb);
}

/**
 * @brief grabs all "long" haplotype-specific nodes, uses tags to find neighbours.
 * @param min_size
 * @param min_ci
 * @param max_ci
 * @return
 */
std::vector<std::vector<std::pair<sgNodeID_t,uint32_t>>> Untangler::find_tag_neighbours(uint32_t min_size, float min_ci, float max_ci) {

    SequenceGraph tsg;
    std::map<sgNodeID_t,sgNodeID_t> general_to_tag;
    sglib::OutputLog()<<"Selecting selected_nodes..."<<std::endl;
    std::vector<std::vector<std::pair<sgNodeID_t,uint32_t>>> neighbours;
    std::vector<std::set<bsg10xTag>> node_tags;
    neighbours.resize(ws.sg.nodes.size());
    node_tags.resize(ws.sg.nodes.size());
    uint64_t total_bp=0;
    auto selected_nodes=ws.select_from_all_nodes(min_size,1000000,20,200000, min_ci, max_ci);
    sglib::OutputLog()<<"Populating node tags..."<<std::endl;
    uint64_t tag_empty=0, tag_10=0, tag_50=0, tag_100=0, tag_1000=0;
    for (auto n:selected_nodes) {
        general_to_tag[n]=tsg.add_node(ws.sg.nodes[n].sequence);
        total_bp+=ws.sg.nodes[n].sequence.size();
        for (auto t:ws.linked_read_mappers[0].get_node_tags(n)) node_tags[n].insert(t);
        if (node_tags[n].empty()) ++tag_empty;
        if (node_tags[n].size()>=10) ++tag_10;
        if (node_tags[n].size()>=50) ++tag_50;
        if (node_tags[n].size()>=100) ++tag_100;
        if (node_tags[n].size()>=1000) ++tag_1000;
    }
    sglib::OutputLog()<<selected_nodes.size()<<" selected totalling "<<total_bp<<"bp "<<std::endl;
    sglib::OutputLog()<<tag_10<<" selected_nodes have 10+ tags"<<std::endl;
    sglib::OutputLog()<<tag_50<<" selected_nodes have 50+ tags"<<std::endl;
    sglib::OutputLog()<<tag_100<<" selected_nodes have 100+ tags"<<std::endl;
    sglib::OutputLog()<<tag_1000<<" selected_nodes have 1000+ tags"<<std::endl;

    sglib::OutputLog()<<"Computing shared tags"<<std::endl;
#pragma omp parallel for shared(neighbours) schedule(static,50)
    for (auto i1=0; i1<selected_nodes.size(); ++i1){
        auto n1=selected_nodes[i1];
        for (auto i2=i1+1;i2<selected_nodes.size();i2++){
            auto n2=selected_nodes[i2];
            uint32_t shared=intersection_size(node_tags[n1],node_tags[n2]);
            if (shared>=10) {
#pragma omp critical
                {
                    if (shared>40) tsg.add_link(general_to_tag[n1],general_to_tag[n2],0);
                    neighbours[n1].emplace_back(n2,shared);
                    neighbours[n2].emplace_back(n1,shared);
                }
            }
        }
    }
    tsg.write_to_gfa("tag_neighbours_nodir.gfa");
    sglib::OutputLog()<<"Sorting shared tags"<<std::endl;
    uint64_t with_neighbours=0;
    for (auto &nn:neighbours){
        if (not nn.empty()) ++with_neighbours;
        std::sort(nn.begin(),nn.end(),[]( const std::pair<sgNodeID_t,uint32_t> &a,
                                          const std::pair<sgNodeID_t,uint32_t> &b) { return a.second>b.second; });
    }
    sglib::OutputLog()<<with_neighbours<<" selected_nodes with neighbours"<<std::endl;
    return neighbours;
}



/**
 * @brief grabs all "long" haplotype-specific nodes, uses tags to find neighbours, uses imbalance to impute direction.
 * @param min_size
 * @param min_ci
 * @param max_ci
 * @return
 */
std::vector<Link>  Untangler::find_tag_neighbours_with_imbalance(uint32_t min_size, float min_ci, float max_ci, float end_perc) {
    //TODO: divide this and make iterative: start by big fish, then add small in-betweens, simplify transitives into subgraphs
    SequenceGraph tsg;
    std::map<sgNodeID_t,sgNodeID_t> general_to_tag,tag_to_general;
    sglib::OutputLog()<<"Selecting selected_nodes..."<<std::endl;
    std::vector<std::vector<std::pair<sgNodeID_t,uint32_t>>> neighbours;
    std::vector<std::set<bsg10xTag>> node_tags;
    neighbours.resize(ws.sg.nodes.size());
    node_tags.resize(ws.sg.nodes.size());
    uint64_t total_bp=0;
    auto selected_nodes=ws.select_from_all_nodes(min_size,1000000,20,200000, min_ci, max_ci);
    for (auto n:selected_nodes) {
        general_to_tag[n]=tsg.add_node(ws.sg.nodes[n].sequence);
        tag_to_general[general_to_tag[n]]=n;
        tag_to_general[-general_to_tag[n]]=-n;
        total_bp+=ws.sg.nodes[n].sequence.size();
        for (auto t:ws.linked_read_mappers[0].get_node_tags(n)) node_tags[n].insert(t);
    }
    sglib::OutputLog()<<selected_nodes.size()<<" selected totalling "<<total_bp<<"bp "<<std::endl;

    sglib::OutputLog()<<"Computing shared tags"<<std::endl;
#pragma omp parallel for shared(tsg) schedule(static,50)
    for (auto i1=0; i1<selected_nodes.size(); ++i1){
        auto n1=selected_nodes[i1];
        for (auto i2=i1+1;i2<selected_nodes.size();i2++){
            auto n2=selected_nodes[i2];
            uint32_t shared=intersection_size(node_tags[n1],node_tags[n2]);
            if (shared>=40) {
//                std::cout<<"Strong connection between "<<n1<<" and "<<n2<<" for tag imbalance"<<std::endl;
                //first create the set of real intersecting tags
                std::set<bsg10xTag> shared_tags;
                std::set_intersection(node_tags[n1].begin(),node_tags[n1].end(),node_tags[n2].begin(),node_tags[n2].end(),std::inserter(shared_tags,shared_tags.end()));
//                std::cout<<"intersection set has "<<shared_tags.size()<<" elements"<<std::endl;
                uint64_t n1_front_in=0,n1_front_total=0,n1_back_in=0,n1_back_total=0;
                uint64_t n2_front_in=0,n2_front_total=0,n2_back_in=0,n2_back_total=0;
                uint64_t n1first30point=ws.sg.nodes[n1].sequence.size()*end_perc;
                uint64_t n1last30point=ws.sg.nodes[n1].sequence.size()*(1-end_perc);
                for (auto rm:ws.linked_read_mappers[0].reads_in_node[n1]){
                    if (rm.first_pos<n1first30point){
                        ++n1_front_total;
                        if (shared_tags.count(ws.linked_read_datastores[0].get_read_tag(rm.read_id))>0) ++n1_front_in;
                    }
                    if (rm.last_pos>n1last30point){
                        ++n1_back_total;
                        if (shared_tags.count(ws.linked_read_datastores[0].get_read_tag(rm.read_id))>0) ++n1_back_in;
                    }
                }
                auto n1f=(100.0*n1_front_in/n1_front_total);
                auto n1b=(100.0*n1_back_in/n1_back_total);
//                std::cout<<" Node "<<n1<<" "<<ws.sg.nodes[n1].sequence.size()<<"bp:  front -> "<<n1f<<" ("<<n1_front_in<<"/"<<n1_front_total<<")"
//                         <<" back -> "<<n1b<<" ("<<n1_back_in<<"/"<<n1_back_total<<")"<<std::endl;
                uint64_t n2first30point=ws.sg.nodes[n2].sequence.size()*end_perc;
                uint64_t n2last30point=ws.sg.nodes[n2].sequence.size()*(1-end_perc);
                for (auto rm:ws.linked_read_mappers[0].reads_in_node[n2]){
                    if (rm.first_pos<n2first30point){
                        ++n2_front_total;
                        if (shared_tags.count(ws.linked_read_datastores[0].get_read_tag(rm.read_id))>0) ++n2_front_in;
                    }
                    if (rm.last_pos>n2last30point){
                        ++n2_back_total;
                        if (shared_tags.count(ws.linked_read_datastores[0].get_read_tag(rm.read_id))>0) ++n2_back_in;
                    }
                }
                auto n2f=(100.0*n2_front_in/n2_front_total);
                auto n2b=(100.0*n2_back_in/n2_back_total);
//                std::cout<<" Node "<<n2<<" "<<ws.sg.nodes[n2].sequence.size()<<"bp:  front -> "<<n2f<<" ("<<n2_front_in<<"/"<<n2_front_total<<")"
//                         <<" back -> "<<n2b<<" ("<<n2_back_in<<"/"<<n2_back_total<<")"<<std::endl;

                if (fabs(2*(n1f-n1b)/(n1f+n1b))>.1 and fabs(2*(n2f-n2b)/(n2f+n2b))>.1) {
//                    std::cout<<"Imbalance solved"<<std::endl;
#pragma omp critical
                    tsg.add_link((n1f > n1b ? general_to_tag[n1] : -general_to_tag[n1]),
                                 (n2f > n2b ? general_to_tag[n2] : -general_to_tag[n2]), 0);
                }
                else {
//                    std::cout<<"Imbalance unsolved"<<std::endl;
                }
//#pragma omp critical
/*                {
                    //if (shared>40) tsg.add_link(general_to_tag[n1],general_to_tag[n2],0);
                    neighbours[n1].emplace_back(n2,shared);
                    neighbours[n2].emplace_back(n1,shared);
                }*/
            }
        }
    }
    tsg.write_to_gfa("tag_neighbours_imbdir.gfa");
    BufferedTagKmerizer btk(ws.linked_read_datastores[0],31,200000,1000);
    uint64_t perf_ts=0;
    for (auto n=1;n<tsg.nodes.size();++n){
        auto b=n;
        auto bfw=tsg.get_fw_links(b);
        auto bbw=tsg.get_bw_links(b);
        if (bbw.size()!=1 or bfw.size()!=1) continue;
        auto a=-bbw[0].dest;
        auto c=bfw[0].dest;
        auto afw=tsg.get_fw_links(a);
        auto cbw=tsg.get_bw_links(c);
        if (afw.size()!=2 or cbw.size()!=2) continue;
        if (afw[0].dest==c or afw[1].dest==c) {
            std::cout << "Found perfect transtition " << a << " " << b << " " << c << std::endl;
            auto A = tag_to_general[a];
            auto B = tag_to_general[b];
            auto C = tag_to_general[c];
            std::cout << "On original ids: " << A << " " << B << " " << C << std::endl;
            std::set<bsg10xTag> tags;
            for (auto t:ws.linked_read_mappers[0].get_node_tags(A)) tags.insert(t);
            for (auto t:ws.linked_read_mappers[0].get_node_tags(B)) tags.insert(t);
            for (auto t:ws.linked_read_mappers[0].get_node_tags(C)) tags.insert(t);
            auto ab = get_all_tag_covered_paths(A, B, tags, btk);
            auto bc = get_all_tag_covered_paths(B, C, tags, btk);
            if (ab.size() == 1) std::cout << "A->B Solved!" << std::endl;
            if (bc.size() == 1) std::cout << "B->C Solved!" << std::endl;
            ++perf_ts;
        }

    }
    std::cout<<perf_ts<<" perfect transitions found"<<std::endl;
    /*sglib::OutputLog()<<"Sorting shared tags"<<std::endl;
    uint64_t with_neighbours=0;
    for (auto &nn:neighbours){
        if (not nn.empty()) ++with_neighbours;
        std::sort(nn.begin(),nn.end(),[]( const std::pair<sgNodeID_t,uint32_t> &a,
                                          const std::pair<sgNodeID_t,uint32_t> &b) { return a.second>b.second; });
    }
    sglib::OutputLog()<<with_neighbours<<" selected_nodes with neighbours"<<std::endl;*/
    return {};
}

uint64_t Untangler::connect_neighbours(uint64_t min_size, float min_ci, float max_ci, int64_t max_distance) {
    //first find all nodes' neighbours
    auto tagneighbours=find_tag_neighbours(min_size,min_ci,max_ci);
    std::vector<std::vector<std::pair<sgNodeID_t,int64_t>>> fndist(ws.sg.nodes.size()),bndist(ws.sg.nodes.size());
    TagWalker tw(ws,{});
    for (auto n=1;n<ws.sg.nodes.size();++n) {
        //explore from a node til hitting a neighbour, check if another selected node is on the way
        std::set<sgNodeID_t> ntn;
        for (auto nd:tagneighbours[n]) ntn.emplace(nd.first);
        if (ntn.empty()) continue;
        fndist[n] = ws.sg.get_distances_to(n, ntn, max_distance);
        bndist[n] = ws.sg.get_distances_to(-n, ntn, max_distance);
    }

    //TODO: only fw or bw neighbour, corresponding, good tag intersection
//    std::cout<<"Finding disconnected tips"<<std::endl;
//    for (auto n=1;n<ws.sg.nodes.size();++n){
//        if (tagneighbours[n].empty()) continue;
//        if (ws.sg.get_fw_links(n).size()==0 or ws.sg.get_bw_links(n).size()==0) {
//            for (auto tn:tagneighbours[n]){
//                if (ws.sg.get_fw_links(tn.first).size()==0 or ws.sg.get_bw_links(tn.first).size()==0) {
//                    std::cout<<"Nodes "<<n<<" and "<<tn.first<<" have disconnected ends and have "<<tn.second<<" shared tags"<<std::endl;
//                }
//            }
//        }
//    }
    //connect_neighbours_trivial(min_size,min_ci,max_ci,max_distance,tagneighbours,bndist,fndist);
    //bndist.resize(ws.sg.nodes.size());
    //fndist.resize(ws.sg.nodes.size());
    //tagneighbours.resize(ws.sg.nodes.size());
    return connect_neighbours_paths_to_same(min_size,min_ci,max_ci,max_distance,tagneighbours,bndist,fndist);
    //std::cout<<"Analysing non-trivial connections"<<std::endl;
    //std::cout


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

void Untangler::connect_neighbours_trivial(uint64_t min_size, float min_ci, float max_ci, int64_t max_distance,
                                           const std::vector<std::vector<std::pair<sgNodeID_t,uint32_t>>> &tagneighbours,
                                           const std::vector<std::vector<std::pair<sgNodeID_t,int64_t>>> & bndist,
                                           const std::vector<std::vector<std::pair<sgNodeID_t,int64_t>>> & fndist) {
    std::cout << "Finding trivial connections" << std::endl;
    std::vector<std::pair<sgNodeID_t, sgNodeID_t>> from_to;
    for (auto n = 1; n < ws.sg.nodes.size(); ++n) {
//todo: check that the other node actually has THIS node as neighbour, validate TAG distances too
        if (fndist[n].size() == 1 and
            (fndist[n][0].first < 0 ? fndist : bndist)[llabs(fndist[n][0].first)].size() == 1) {
//std::cout<<"Nodes "<<n<<" and "<<fndist[n][0].first<<" may be trivial neighbours"<<std::endl;
            if (llabs(fndist[n][0].first) > n) from_to.push_back({n, fndist[n][0].first});
        }
        if (bndist[n].size() == 1 and
            (bndist[n][0].first < 0 ? fndist : bndist)[llabs(bndist[n][0].first)].size() == 1) {
//std::cout<<"Nodes "<<-n<<" and "<<bndist[n][0].first<<" may be trivial neighbours"<<std::endl;
            if (llabs(bndist[n][0].first) > n) from_to.push_back({-n, bndist[n][0].first});
        }
//should we check also that there is only one single connection out and in the nodes?

    }
    std::cout << from_to.size() << " trivial connections to inflate" << std::endl;
    uint64_t done = 0;
    for (auto ft:from_to) {
//Create all possible paths between from_to (sg function):
        auto allpaths = ws.sg.find_all_paths_between(ft.first, ft.second, 50000);
//  is there only one? does it contain no selected node? other conditions? Optional: check the path is covered by kmers from tags on both ends of it.
        if (allpaths.size() != 1) continue;
        for (auto &n: allpaths[0].nodes) if (not tagneighbours[llabs(n)].empty()) continue;
        if (allpaths[0].nodes.empty()) continue; //XXX:these should actually be connected directly?
//  create a single node with the sequence of the path, migrate connections from the nodes (sg function)
        auto new_node = ws.sg.add_node(Node(allpaths[0].get_sequence()));
        auto old_fw = ws.sg.get_fw_links(ft.first);
        for (auto &l:old_fw) {
            if (l.dest == allpaths[0].nodes.front()) ws.sg.add_link(-ft.first, new_node, l.dist);
            ws.sg.remove_link(l.source, l.dest);
        }
        auto old_bw = ws.sg.get_bw_links(ft.second);
        for (auto &l:old_bw) {
            if (l.dest == -allpaths[0].nodes.back()) ws.sg.add_link(ft.second, -new_node, l.dist);
            ws.sg.remove_link(l.source, l.dest);
        }
//std::cout<<"Replaced connection between "<<ft.first<<" and "<<ft.second<<" through nodes ";
//for (auto &n:allpaths[0].nodes) std::cout<<n<<" ";
//std::cout<<" by a single node "<<new_node<<std::endl;
        ++done;


    }
    std::cout << done << " connections inflated" << std::endl;
}

uint64_t Untangler::connect_neighbours_paths_to_same(uint64_t min_size, float min_ci, float max_ci, int64_t max_distance,
                                           const std::vector<std::vector<std::pair<sgNodeID_t,uint32_t>>> &tagneighbours,
                                           const std::vector<std::vector<std::pair<sgNodeID_t,int64_t>>> & bndist,
                                           const std::vector<std::vector<std::pair<sgNodeID_t,int64_t>>> & fndist) {
    std::cout << "Finding multi-path connections to same neighbour" << std::endl;
    std::vector<std::pair<sgNodeID_t, sgNodeID_t>> from_to;
    std::vector<sgNodeID_t> single_fw(ws.sg.nodes.size(),0),single_bw(ws.sg.nodes.size(),0);
    //first mark outputs going into a single neighbour
    for (auto n = 1; n < ws.sg.nodes.size(); ++n) {
        if (fndist[n].size()>=1) {
            auto fdest = fndist[n][0].first;
            bool fok = true;
            for (auto fd:fndist[n]) if (fdest != fd.first) fok = false;
            if (fok) single_fw[n] = fdest;
        }
        if (bndist[n].size()>=1) {
            auto bdest = bndist[n][0].first;
            bool bok = true;
            for (auto bd:bndist[n]) if (bdest != bd.first) bok = false;
            if (bok) single_bw[n] = bdest;
        }
    }
    //now look for happy matches where both have a single connection and it is the same one!
    for (auto n = 1; n < ws.sg.nodes.size(); ++n) {
        if (single_fw[n]>0 and single_fw[n]<n and single_bw[single_fw[n]]==-n) from_to.push_back({n,single_fw[n]});
        if (single_fw[n]<0 and -single_fw[n]<n and single_fw[-single_fw[n]]==-n) from_to.push_back({n,single_fw[n]});
        if (single_bw[n]>0 and single_bw[n]<n and single_bw[single_bw[n]]==n) from_to.push_back({-n,single_bw[n]});
        //std::cout<<"."<<std::flush;
        if (single_bw[n]<0 and -single_bw[n]<n and single_fw[-single_bw[n]]==n) from_to.push_back({-n,single_bw[n]});
        //std::cout<<"."<<std::endl;
    }
    std::cout << from_to.size() << " multi-path connections to same neighbour to inflate" << std::endl;

    uint64_t done = 0;
    std::vector<std::pair<std::pair<sgNodeID_t,sgNodeID_t>,SequenceGraphPath>> final_sols;
#pragma omp parallel
    {
        BufferedTagKmerizer btk(ws.linked_read_datastores[0],31,200000,1000);
#pragma omp for schedule(dynamic,10)
        for (auto idx = 0; idx < from_to.size(); ++idx) {
            auto &ft = from_to[idx];
//Create all possible paths between from_to (sg function):
            auto allpaths = ws.sg.find_all_paths_between(ft.first, ft.second, max_distance);

            bool complex = false;
            for (auto p:allpaths) {
                for (auto &n: p.nodes) {
                    if (not tagneighbours[llabs(n)].empty()) {
                        complex = true;
                        //std::cout<<"A path node has other neighbours!"<<std::endl;
                    }
                }
                if (p.nodes.empty()) {
                    complex = true;
                    //std::cout<<"Path is empty!"<<std::endl;
                }
            }
            if (complex or allpaths.empty()) continue;
            std::cout << "Evaluating between " << allpaths.size() << " different paths to join " << ft.first << " and "
                      << ft.second << std::endl;
            auto ntags = ws.linked_read_mappers[0].get_node_tags(ft.first);
            for (auto t:ws.linked_read_mappers[0].get_node_tags(ft.second)) ntags.insert(t);
            auto neightagkmers = btk.get_tags_kmers(6, ntags);
            std::vector<SequenceGraphPath> sol;
            StringKMerFactory skf(31);
            for (auto p:allpaths) {
                auto pseq = p.get_sequence();
                //TODO: this can be doen not with a factory but with a class that already counts how many are in the set
                std::vector<uint64_t> seqkmers;
                skf.create_kmers(pseq,seqkmers);
                uint64_t covered = 0, uncovered = 0;
                for (auto x:seqkmers) {
                    if (neightagkmers.count(x) > 0) ++covered;
                    else ++uncovered;
                }
                std::cout << "  path evaluation: " << covered << " covered, " << uncovered << " uncovered" << std::endl;
                if (uncovered == 0) sol.push_back(p);
            }
            if (sol.size() == 1) {
                ++done;
                std::cout << "JOIN SOLVED!" << std::endl;
#pragma omp critical
                final_sols.emplace_back(ft, sol[0]);
            }


        }
    }
    std::cout << done << " joins solved (connections to inflate)" << std::endl;
    for (auto &sol:final_sols){
        dettach_path_as_new_node(sol.first.first,sol.first.second,sol.second);
    }
    return final_sols.size();
}

void Untangler::dettach_path_as_new_node(sgNodeID_t from, sgNodeID_t to, SequenceGraphPath path) {
    auto new_node = ws.sg.add_node(Node(path.get_sequence()));
    auto old_fw = ws.sg.get_fw_links(from);
    for (auto &l:old_fw) {
        if (l.dest == path.nodes.front()) ws.sg.add_link(-from, new_node, l.dist);
        ws.sg.remove_link(l.source, l.dest);
    }
    auto old_bw = ws.sg.get_bw_links(to);
    for (auto &l:old_bw) {
        if (l.dest == -path.nodes.back()) ws.sg.add_link(to, -new_node, l.dist);
        ws.sg.remove_link(l.source, l.dest);
    }
}

std::vector<SequenceGraphPath> Untangler::get_all_tag_covered_paths(sgNodeID_t from, sgNodeID_t to,
                                                                    std::set<bsg10xTag> &tags,
                                                                    BufferedTagKmerizer &btk) {
    int64_t max_distance=50000;
    auto allpaths = ws.sg.find_all_paths_between(from, to, max_distance);


    if (allpaths.empty()) return {};
    std::cout << "Evaluating between " << allpaths.size() << " different paths to join " << from << " and "
              << to << std::endl;
    auto neightagkmers = btk.get_tags_kmers(6, tags);
    std::vector<SequenceGraphPath> sol;
    StringKMerFactory skf(31);
    for (auto p:allpaths) {
        auto pseq = p.get_sequence();
        //TODO: this can be doen not with a factory but with a class that already counts how many are in the set
        std::vector<uint64_t> seqkmers;
        skf.create_kmers(pseq, seqkmers);
        uint64_t covered = 0, uncovered = 0;
        for (auto x:seqkmers) {
            if (neightagkmers.count(x) > 0) ++covered;
            else ++uncovered;
        }
        std::cout << "  path evaluation: " << covered << " covered, " << uncovered << " uncovered" << std::endl;
        if (uncovered == 0) sol.push_back(p);
    }
    return sol;
}

/**
 * @brief returns the percentage of reads in both ends of node covered by tags in tags
 * @param node
 * @param tags
 * @param end_perc
 * @param end_size
 * @return
 */
std::pair<float,float> Untangler::tag_read_percentage_at_ends(sgNodeID_t node, std::set<bsg10xTag> tags, float end_perc,
                                                              uint32_t end_size) {
    auto n1=llabs(node);
    if (end_size==0) end_size=ws.sg.nodes[n1].sequence.size()*end_perc;
    auto first_chunk=end_size;
    auto last_chunk=ws.sg.nodes[n1].sequence.size()-end_size;
    uint64_t n1_front_in=0,n1_front_total=0,n1_back_in=0,n1_back_total=0;
    for (auto rm:ws.linked_read_mappers[0].reads_in_node[n1]){
        if (rm.first_pos<first_chunk){
            ++n1_front_total;
            if (tags.count(ws.linked_read_datastores[0].get_read_tag(rm.read_id))>0) ++n1_front_in;
        }
        if (rm.last_pos>last_chunk){
            ++n1_back_total;
            if (tags.count(ws.linked_read_datastores[0].get_read_tag(rm.read_id))>0) ++n1_back_in;
        }
    }
    auto n1f=(100.0*n1_front_in/n1_front_total);
    auto n1b=(100.0*n1_back_in/n1_back_total);
    return {n1f,n1b};
}

/**
 * @brief, creates backbones by joining tag-neighbours with tag imbalance. starts by larger nodes
 * - Every selected node must appear in only one backbone.
 * - Nodes can be removed from selection.
 * - When incorporating smaller nodes, the selection for large nodes remains unchanged.
 * @param min_ci
 * @param max_ci
 * @param end_perc
 * @return
 */
std::vector<Backbone> Untangler::create_backbones(uint64_t min_size, float min_ci, float max_ci, float end_perc, int min_shared_tags) {
    std::vector<Backbone> backbones;
    auto selected_nodes=ws.select_from_all_nodes(min_size,1000000,20,200000, min_ci, max_ci);
    std::vector<std::pair<sgNodeID_t , sgNodeID_t>> linked_nodes;
    linked_nodes.reserve(10*selected_nodes.size());
    std::vector<std::set<bsg10xTag >> selected_nodes_tags;
    selected_nodes_tags.reserve(selected_nodes.size());
    std::cout<<"Computing directional connections among "<<selected_nodes.size()<<" nodes"<<std::endl;
    for (auto n:selected_nodes) selected_nodes_tags.push_back(ws.linked_read_mappers[0].get_node_tags(n));
#pragma omp parallel for schedule (static,50)
    for (auto i1=0;i1<selected_nodes.size();++i1){
        auto n1=selected_nodes[i1];
        for (auto i2=i1+1;i2<selected_nodes.size();++i2){
            auto n2=selected_nodes[i2];
            std::set<bsg10xTag> shared_tags;
            std::set_intersection(selected_nodes_tags[i1].begin(),selected_nodes_tags[i1].end(),
                                  selected_nodes_tags[i2].begin(),selected_nodes_tags[i2].end(),
                                  std::inserter(shared_tags,shared_tags.end()));
            if (shared_tags.size()>min_shared_tags) {
                auto n1fb=tag_read_percentage_at_ends(n1,shared_tags, end_perc);
                auto n1f=n1fb.first; auto n1b=n1fb.second;
                auto n2fb=tag_read_percentage_at_ends(n2,shared_tags, end_perc);
                auto n2f=n2fb.first; auto n2b=n2fb.second;

                if (fabs(2 * (n1f - n1b) / (n1f + n1b)) > .1 and fabs(2 * (n2f - n2b) / (n2f + n2b)) > .1) {
                    //std::cout<<"Imbalance solved with "<<shared_tags.size()<<" tags "<<n1<< "( "<<n1f<<", "<<n1b<<" ) "<<n2<< "( "<<n2f<<", "<<n2b<<" )"<<std::endl;
#pragma omp critical
                    {
                        linked_nodes.emplace_back((n1f > n1b ? -i1 : i1), (n2f > n2b ? -i2 : i2));
                    }

                }
            }

        }
    }
    std::cout<<"Total connections detected: "<<linked_nodes.size()<<std::endl;

    //transform links into fw and bw connections (unelegant!)
    std::vector<std::vector<sgNodeID_t>> fw_neighbours(selected_nodes.size());
    std::vector<std::vector<sgNodeID_t>> bw_neighbours(selected_nodes.size());
    for (auto l:linked_nodes) {
        if (l.first>0) bw_neighbours[l.first].push_back(l.second);
        else fw_neighbours[-l.first].push_back(l.second);
        if (l.second>0) bw_neighbours[l.second].push_back(l.first);
        else fw_neighbours[-l.second].push_back(l.first);
    }

    std::ofstream of_node_connections("backbone_node_connections.csv"), of_links("backbone_links.csv");
    for (auto l:linked_nodes) of_links<<(l.first>0 ? selected_nodes[l.first]:-selected_nodes[-l.first])<<", "
                                      <<(l.second>0 ? selected_nodes[l.second]:-selected_nodes[-l.second])<<std::endl;
    // check nodes with only one connection
    uint64_t none=0, one_side=0, one_in_one_out=0, multiple=0, two=0;
    for (auto i1=0;i1<selected_nodes.size();++i1) {
        auto n1 = selected_nodes[i1];
        if (fw_neighbours[i1].size()==0 and bw_neighbours[i1].size()==0) ++none;
        else if (fw_neighbours[i1].size()==0 or bw_neighbours[i1].size()==0) ++one_side;
        else if (fw_neighbours[i1].size()==1 and bw_neighbours[i1].size()==1) ++one_in_one_out;
        else if (fw_neighbours[i1].size()==2 or bw_neighbours[i1].size()==2) ++two;
        else ++multiple;

        of_node_connections<<n1;
        for (auto i2:bw_neighbours[i1]) of_node_connections<<", "<<(i2>0 ? selected_nodes[i2]:-selected_nodes[-i2]);
        of_node_connections<<std::endl;
        of_node_connections<<-n1;
        for (auto i2:fw_neighbours[i1]) of_node_connections<<", "<<(i2>0 ? selected_nodes[i2]:-selected_nodes[-i2]);
        of_node_connections<<std::endl;
    }
    std::cout<<none<<" nodes no in or out connections"<<std::endl;
    std::cout<<one_side<<" nodes have one-side-only connections"<<std::endl;
    std::cout<<one_in_one_out<<" nodes have 1 in and 1 out connections"<<std::endl;
    std::cout<<two<<" nodes have 2 connections on at least one side"<<std::endl;
    std::cout<<multiple<<" nodes have multiple connections, with at least 1 per side"<<std::endl;


    //check and remove simple linear transitive connections

    //check nodes with multiple incoherent transtitions (repeats)
/*
    std::vector<bool> used(selected_nodes.size(),false);
    for (auto i1=0;i1<selected_nodes.size();++i1) {
        if (selected_nodes[i1]) continue;
        auto n1 = selected_nodes[i1];
        backbones.emplace_back(ws,*this);
        used[n1]=true;
        std::vector<sgNodeID_t> last_added={i1};
        backbones.back().nodes.push_back(n1);
        while (not last_added.empty()){
            std::vector<sgNodeID_t> new_last_added;
            for (auto x:last_added){
                for (auto y:node_neighbours[x]){
                    if (not used[y]){
                        used[y]=true;
                        backbones.back().nodes.push_back(selected_nodes[y]);
                        new_last_added.push_back(y);
                    }
                }
            }
            last_added=new_last_added;
        }
        std::cout<<"new backbone: ";
        for (auto n:backbones.back().nodes) std::cout<<"seq"<<n<<", ";
        std::cout<<std::endl;
    }
*/
    return backbones;

}

void Untangler::unroll_simple_loops() {
    std::vector<bool> used(ws.sg.nodes.size(),false);
    for (auto i=1;i<ws.sg.nodes.size();++i){
        if (ws.sg.nodes[i].status==sgNodeDeleted) continue;
        if (used[i]) continue;
        auto fwl=ws.sg.get_fw_links(i);
        auto bwl=ws.sg.get_bw_links(i);
        for (auto f:fwl) {
            for (auto b:bwl) {
                if (f.dest == -b.dest) {
                    used[i]=true;
                    used[llabs(f.dest)]=true;
                    //std::cout<< "simple loop detected between "
                    //         <<i<< " ("<<ws.sg.nodes[i].sequence.size()<<"bp, kci="<<ws.kci.compute_compression_for_node(i,1)<<") and "
                    //         << f.dest <<" ("<<ws.sg.nodes[llabs(f.dest)].sequence.size()<<"bp, kci="<<ws.kci.compute_compression_for_node(llabs(f.dest),1)<<")"<<std::endl;
                    std::cout<<"seq"<<i<< ", "<<ws.sg.nodes[i].sequence.size()<<", "<<ws.kci.compute_compression_for_node(i,1)<<", seq"
                            << llabs(f.dest) <<", "<<ws.sg.nodes[llabs(f.dest)].sequence.size()<<", "<<ws.kci.compute_compression_for_node(llabs(f.dest),1)<<std::endl;
                }
            }
        }

    }
}

void PairedReadLinker::generate_links_size_ci( uint32_t min_size, float min_ci, float max_ci,int min_reads) {
    std::vector<bool> to_link(ws.sg.nodes.size());
    for (auto &n:ws.select_from_all_nodes(min_size,1000000000,0,1000000000,min_ci,max_ci)) to_link[n]=true;
    generate_links(to_link,min_reads);
}
void PairedReadLinker::generate_links_hspnp(int min_reads) {
    std::vector<bool> to_link(ws.sg.nodes.size());
    uint64_t count=0,bp=0;
    for (auto &n:u.get_all_HSPNPs()) {
        to_link[llabs(n.first)]=true;
        to_link[llabs(n.second)]=true;
        ++count;
        bp+=ws.sg.nodes[llabs(n.first)].sequence.size()+ws.sg.nodes[llabs(n.second)].sequence.size();
    }
    sglib::OutputLog()<<count<<" nodes in HSPNP selected for linkage totalling "<<bp<<std::endl;
    generate_links(to_link,min_reads);
}
void PairedReadLinker::generate_links( const std::vector<bool> &to_link,int min_reads) {

    sglib::OutputLog()<<"filling orientation indexes"<<std::endl;
    uint64_t revc=0,dirc=0,false_rev=0,false_dir=0,true_rev=0,true_dir=0;
    std::vector<std::vector<bool>> orientation;
    for (auto &pm:ws.paired_read_mappers){
        orientation.emplace_back();
        orientation.back().resize(pm.read_to_node.size());
        for (auto n=1;n<ws.sg.nodes.size();++n)
            for (auto &rm:pm.reads_in_node[n]) {
            orientation.back()[rm.read_id]=rm.rev;
            if (rm.first_pos<rm.last_pos){if (rm.rev) ++false_rev; else ++true_rev;};
            if (rm.first_pos>rm.last_pos ){if (!rm.rev) ++false_dir; else ++true_dir;};
            if (rm.rev) revc++;
            else dirc++;
        }
    }
    std::ofstream lof("paired_links.txt");
    sglib::OutputLog()<<"FW: "<<dirc<<" ( "<<true_dir<<" - "<< false_dir<<" )"<<std::endl;
    sglib::OutputLog()<<"BW: "<<revc<<" ( "<<true_rev<<" - "<< false_rev<<" )"<<std::endl;
    std::map<std::pair<sgNodeID_t, sgNodeID_t>, uint64_t> lv;
    sglib::OutputLog()<<"collecting link votes across all paired libraries"<<std::endl;
    //use all libraries collect votes on each link
    auto rmi=0;
    for (auto &pm:ws.paired_read_mappers) {
        for (auto i = 1; i < pm.read_to_node.size(); i += 2) {
            sgNodeID_t n1 = pm.read_to_node[i];
            sgNodeID_t n2 = pm.read_to_node[i + 1];
            if (n1 == 0 or n2 == 0 or n1 == n2 or !to_link[n1] or !to_link[n2]) continue;
            if (orientation[rmi][i]) n1=-n1;
            if (orientation[rmi][i+1]) n2=-n2;
            if (llabs(n1) > llabs(n2)) std::swap(n1,n2);
            ++lv[std::make_pair(n1, n2)];
        }
        ++rmi;
    }

    sglib::OutputLog()<<"adding links with "<<min_reads<<" votes"<<std::endl;
    //std::vector<std::vector<std::pair<sgNodeID_t ,uint64_t>>> nodelinks(ws.sg.nodes.size());
    for (auto l:lv) {
        if (l.second>=min_reads){
            //todo: size, appropriate linkage handling, etc
            //todo: check alternative signs for same linkage
            auto s=l.first.first;
            auto d=l.first.second;
            auto v1=std::make_pair(-s,d);
            auto v2=std::make_pair(-s,-d);
            auto v3=std::make_pair(s,-d);
            if (lv.count(v1) and lv[v1]>5*l.second) continue;
            if (lv.count(v2) and lv[v2]>5*l.second) continue;
            if (lv.count(v3) and lv[v3]>5*l.second) continue;
            add_link(l.first.first,l.first.second,0);
            //lof<<l.first.first<<" "<<l.first.second<<" "<<l.second<<std::endl;
        }
    }
    sglib::OutputLog()<<"dumping links"<<std::endl;
    for (auto n=1;n<ws.sg.nodes.size();++n) {
        auto fwl=get_fw_links(n);
        auto bwl=get_bw_links(n);
        if (!fwl.empty()) {
            lof<<n<<" FW: ";
            for (auto &l:fwl) lof<<" "<<l.dest;
            lof<<std::endl<<"          seq"<<n;
            for (auto &l:fwl) lof<<", seq"<<llabs(l.dest);
            lof<<std::endl;
        }
        if (!bwl.empty()) {
            lof<<n<<" BW: ";
            for (auto &l:bwl) lof<<" "<<l.dest;
            lof<<std::endl<<"          seq"<<n;
            for (auto &l:bwl) lof<<", seq"<<llabs(l.dest);
            lof<<std::endl;
        }
    }
    //TODO: remove printing
//    lof<<n<<": ";
//    for (auto l:flv) lof<<" "<<l.first<<"("<<l.second<<")";
//    lof<<std::endl;
    //filter by min reads
    //infer size?
    //TODO:hardcoded LMP orientation


}

void PairedReadLinker::add_link(sgNodeID_t source, sgNodeID_t dest, int32_t d) {
    Link l(source,dest,d);
    links[(source > 0 ? source : -source)].emplace_back(l);
    std::swap(l.source,l.dest);
    links[(dest > 0 ? dest : -dest)].emplace_back(l);
}

void PairedReadLinker::remove_link(sgNodeID_t source, sgNodeID_t dest) {
    auto & slinks = links[(source > 0 ? source : -source)];
    slinks.erase(std::remove(slinks.begin(), slinks.end(), Link(source,dest,0)), slinks.end());
    auto & dlinks = links[(dest > 0 ? dest : -dest)];
    dlinks.erase(std::remove(dlinks.begin(), dlinks.end(), Link(dest,source,0)), dlinks.end());

}

std::vector<Link> PairedReadLinker::get_fw_links( sgNodeID_t n){
    std::vector<Link> r;
    for (auto &l:links[(n>0 ? n : -n)]) if (l.source==-n) r.emplace_back(l);
    return r;
}

//returns a list of all fw nodes up to radius jumps away.
std::set<sgNodeID_t> PairedReadLinker::fw_reached_nodes(sgNodeID_t n, int radius) {
    std::set<sgNodeID_t> reached,last={n};
    for (auto i=0;i<radius;++i) {
        std::set<sgNodeID_t> new_last;
        for (auto l:last){
            for (auto fwl:get_fw_links(l)){
                new_last.insert(fwl.dest);
                reached.insert(fwl.dest);
            }
        }
        std::swap(last,new_last);
    }
    return reached;
}

void PairedReadLinker::remove_transitive_links(int radius) {
    std::cout<<"removing transitive connections..."<<std::endl;

    //TODO: check mutually transitive connections (A->B, B->C, A->C, C->B)
    for (auto fn=1;fn<ws.sg.nodes.size();++fn) {
        for (auto n:{fn,-fn}) {
            //create a list of fw nodes reached through each fw connection in up to radius steps

            auto fwl = get_fw_links(n);
            std::vector<sgNodeID_t> neighbours;
            std::vector<std::set<sgNodeID_t>> reaches;
            for (auto fl:fwl){
                neighbours.push_back(fl.dest);
                reaches.push_back(fw_reached_nodes(fl.dest,radius));
            }
            if (neighbours.empty()) continue;
            std::set<sgNodeID_t> indirect;
            for (auto bi=0;bi<neighbours.size();++bi){
                auto b=neighbours[bi];
                for (auto ci=0;ci<neighbours.size();++ci){
                    if (bi==ci) continue;
                    auto c=neighbours[ci];
                    //if C is reached through B but B is not through C, then remove connection to C.
                    if (reaches[bi].count(c)>0 and reaches[ci].count(b)==0){
                        //std::cout << "Transitive connection detected! " << n << " -> " << b << " -> "
                        //          << c << std::endl;
                        //std::cout << "  seq" << llabs(n) << ", seq" << b << ", seq" << llabs(c) << std::endl;
                        indirect.insert(c);
                    }
                }

            }
            std::cout<<"Connections for node " <<labs(n) << (n>0 ? " FW":" BW")<<":"<<std::endl;
            for (auto l:neighbours) {
                std::cout<<l;
                if (indirect.count(l)>0) std::cout<<" indirect"<<std::endl;
                else std::cout<<" DIRECT!"<<std::endl;
            }
        }

    }
    std::cout<<"removing transitive connections DONE"<<std::endl;
}

std::vector<std::vector<sgNodeID_t>> PairedReadLinker::find_local_problems(uint64_t long_node_size) {
    std::vector<std::vector<sgNodeID_t>> local_problem;
    std::vector<bool> used_fw(ws.sg.nodes.size());
    std::vector<bool> used_bw(ws.sg.nodes.size());
    for (auto n=1;n<ws.sg.nodes.size();++n) {
        //start from a long node going fw (and later bw)
        if (ws.sg.nodes[n].sequence.size()<long_node_size) continue;
        for (auto sn:{n,-n}) {
            if (sn==n and used_fw[n]) continue;
            if (sn==-n and used_bw[n]) continue;
            //get all bidirectional neighbours from all other nodes, dont go through long nodes (frontiers)
            //std::cout<<"checking large node "<<n<<std::endl;
            std::set<sgNodeID_t> internal_nodes;
            std::set<sgNodeID_t> frontiers={sn};
            bool mod=true;
            while (mod) {
                //std::cout<<" new expansion round "<<std::endl;
                //std::cout<<"  frontiers: ";
                //for (auto f:frontiers)std::cout<<" "<<f;
                //std::cout<<std::endl;
                //std::cout<<"  internals: ";
                //for (auto f:internal_nodes)std::cout<<" "<<f;
                //std::cout<<std::endl;
                mod=false;
                //grab all neighbours of internal nodes and all fw_neighbours of frontiers till no new nodes or too many
                auto fwn=frontiers;
                for (auto in:internal_nodes) {
                    fwn.insert(in);
                    fwn.insert(-in);
                }
                for (auto f:fwn) {
                    //std::cout<<"  checking fw links of "<<f<<std::endl;
                    for (auto fwl:get_fw_links(f)) {
                        auto en=fwl.dest;
                        //std::cout<<"    found connection to "<<en<<std::endl;
                        if (ws.sg.nodes[llabs(en)].sequence.size()<long_node_size){
                            //std::cout<<"     "<<en<<" is short"<<std::endl;
                            if (internal_nodes.count(en)==0 and internal_nodes.count(-en)==0){
                                mod=true;
                                internal_nodes.insert(en);
                                //std::cout<<"     "<<en<<" added to internal nodes"<<std::endl;
                            }
                        }
                        else {
                            //std::cout<<"     "<<en<<" is long"<<std::endl;
                            if (frontiers.count(-en)==0){
                                mod=true;
                                frontiers.insert(-en);
                                //std::cout<<"     "<<-en<<" added to frontiers"<<std::endl;
                            }
                        }
                    }
                }
                if (frontiers.size()+internal_nodes.size()>100) break;
            }
            //add result to local_problem
            if (frontiers.size()+internal_nodes.size()>100) continue;
            local_problem.push_back({});
            for (auto f:frontiers) {
                local_problem.back().push_back(f);
                if (f>0) used_fw[f]=true;
                else used_bw[-f]=true;
            }
            for (auto f:internal_nodes) local_problem.back().push_back(f);
        }
    }
    return local_problem;
}

std::vector<std::vector<sgNodeID_t>> PairedReadLinker::solve_local_problem(std::vector<sgNodeID_t> connected_nodes) {
    PairedReadLinker local_linker(ws, u, *this, connected_nodes); //create a local linker with only the selected node's connections
    local_linker.remove_transitive_links(10);
    //if (local_linker.is_fully_solved()) return local_linker.get_all_lines();
    return {};

}