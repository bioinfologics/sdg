//
// Created by Bernardo Clavijo (EI) on 26/02/2018.
//

#include "Untangler.hpp"
#include "TagWalker.hpp"


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
    std::atomic<uint64_t > processing(0);
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