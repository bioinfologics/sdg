//
// Created by Bernardo Clavijo (EI) on 28/05/2018.
//

#include "LinkageUntangler.hpp"

void LinkageUntangler::clear_node_selection() {
    selected_nodes.clear();
    selected_nodes.resize(ws.sg.nodes.size());
    frontier_nodes.clear();
    frontier_nodes.resize(ws.sg.nodes.size());
}

void LinkageUntangler::report_node_selection() {
    uint64_t total_bp=0,total_count=0,selected_bp=0,selected_count=0;
    for (auto n=1;n<ws.sg.nodes.size();++n) {
        if (ws.sg.nodes[n].status == sgNodeDeleted) continue;
        total_bp+=ws.sg.nodes[n].sequence.size();
        ++total_count;
        if (selected_nodes[n]) {
            selected_bp += ws.sg.nodes[n].sequence.size();
            ++selected_count;
        }
    }
        sglib::OutputLog()<< "Current selection: "<<selected_count<<" / "<<total_count<<" nodes  with  "<<selected_bp<<" / "<<total_bp<<" bp"<<std::endl;

}

void LinkageUntangler::select_nodes_by_size_and_ci( uint64_t min_size, float min_ci, float max_ci) {
    std::vector<sgNodeID_t> nodes;
    sglib::OutputLog()<<"LU selecting nodes by size and ci: size >= " << min_size << " bp  |  " << min_ci << "<= CI <=" << max_ci <<std::endl;
#pragma omp parallel
    {
#pragma omp for schedule(static, 100)
        for (auto n=1;n<ws.sg.nodes.size();++n) {
            if (ws.sg.nodes[n].status==sgNodeDeleted) continue;
            if (ws.sg.nodes[n].sequence.size() < min_size) continue;
            auto ci = ws.kci.compute_compression_for_node(n, 1);
            if (std::isnan(ci) or ci < min_ci or ci > max_ci) continue;
            #pragma omp critical(collect_selected_nodes)
            selected_nodes[n]=true;
        }

    }
}

void LinkageUntangler::select_nodes_by_HSPNPs(uint64_t min_size, float min_ci, float max_ci) {
    //first, get all perfectly parallel nodes (TODO:generalise more)
    //std::ofstream hl("hspnp_list.txt");
    std::set<std::pair<sgNodeID_t, sgNodeID_t >> pnps, hspnps;
#pragma omp parallel for schedule(static, 100)
    for (sgNodeID_t n = 1; n < ws.sg.nodes.size(); ++n) {
        if (ws.sg.nodes[n].status == sgNodeDeleted) continue;
        if (ws.sg.nodes[n].sequence.size() < min_size) continue;
        //FW check
        auto fwl = ws.sg.get_fw_links(n);
        if (fwl.size() != 1) continue;
        auto post = fwl[0].dest;
        auto post_bwl = ws.sg.get_bw_links(post);
        if (post_bwl.size() != 2) continue;
        if (llabs(post_bwl[0].dest)==llabs(post_bwl[1].dest))continue;
        //BW check
        auto bwl = ws.sg.get_bw_links(n);
        if (bwl.size() != 1) continue;
        auto prev = bwl[0].dest;
        auto prev_fwl = ws.sg.get_bw_links(prev);
        if (prev_fwl.size() != 2) continue;

        if ((prev_fwl[0].dest == -post_bwl[0].dest and prev_fwl[1].dest == -post_bwl[1].dest)
            or (prev_fwl[1].dest == -post_bwl[0].dest and prev_fwl[0].dest == -post_bwl[1].dest)) {
            sgNodeID_t m;
            if (llabs(prev_fwl[0].dest) != n and llabs(prev_fwl[1].dest) != n) std::cout<<"Error! cant find N in prev!"<<std::endl;
            if (llabs(prev_fwl[0].dest) == n) m = llabs(prev_fwl[1].dest);
            else m = llabs(prev_fwl[0].dest);
#pragma omp critical(inserting_pnps)
            {
                if (n < m) pnps.insert(std::make_pair(n, m));
                else pnps.insert(std::make_pair(m, n));
            }
            //Now evaluate coverage of the branches
            auto c1 = ws.kci.compute_compression_for_node(n, 1);
            if (std::isnan(c1) or c1<min_ci or c1>max_ci) continue;
            auto c2 = ws.kci.compute_compression_for_node(m, 1);
            if (std::isnan(c2) or c2<min_ci or c2>max_ci) continue;
#pragma omp critical(inserting_hspnps)
            {
                //hl<<(n<m ? n:m)<<" "<<(n<m ? m:n)<<std::endl;
                if (n < m) hspnps.insert(std::make_pair(n, m));
                else hspnps.insert(std::make_pair(m, n));
            }
        }
    }

    sglib::OutputLog() << "Selecting HSPNPs: " << pnps.size() << " pairs passed topology, " << hspnps.size()
                       << " passed CI" << std::endl;
    for (auto p:hspnps) {
        selected_nodes[p.first] = true;
        selected_nodes[p.second] = true;
    }

}

LinkageDiGraph LinkageUntangler::make_topology_linkage(int radius) {
    LinkageDiGraph ldg(ws.sg);
    for (auto m=1;m<ws.sg.nodes.size();++m) {
        if (!selected_nodes[m]) continue;
        for (auto n:{m,-m}) {
            std::set<sgNodeID_t> reached, last = {n};
            for (auto i = 0; i < radius; ++i) {
                std::set<sgNodeID_t> new_last;
                for (auto l:last) {
                    for (auto fwl:ws.sg.get_fw_links(l)) {
                        if (selected_nodes[llabs(fwl.dest)]) {
                            ldg.add_link(-n, fwl.dest, 0);
                        } else {
                            new_last.insert(fwl.dest);
                        }

                    }
                }
                std::swap(last, new_last);
            }
        }
    }
    return ldg;
}

LinkageDiGraph LinkageUntangler::make_paired_linkage(int min_reads) {
    LinkageDiGraph ldg(ws.sg);
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
            if (n1 == 0 or n2 == 0 or n1 == n2 or !selected_nodes[n1] or !selected_nodes[n2] ) continue;
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
            ldg.add_link(l.first.first,l.first.second,0);
            //lof<<l.first.first<<" "<<l.first.second<<" "<<l.second<<std::endl;
        }
    }
    return ldg;
}