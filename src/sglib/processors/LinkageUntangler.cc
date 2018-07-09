//
// Created by Bernardo Clavijo (EI) on 28/05/2018.
//

#include "LinkageUntangler.hpp"
#include "GraphEditor.hpp"
#include "GraphMaker.hpp"

class KmerMapCreator : public  KMerFactory {
public:
    explicit KmerMapCreator(uint8_t k) : KMerFactory(k){}
    inline void create_all_kmers(const char * seq, std::unordered_map<uint64_t,uint32_t> &mers){
        // TODO: Adjust for when K is larger than what fits in uint64_t!
        last_unknown=0;
        fkmer=0;
        rkmer=0;
        auto s=seq;
        while (*s!='\0' and *s!='\n') {
            //fkmer: grows from the right (LSB)
            //rkmer: grows from the left (MSB)
            fillKBuf(*s, 0, fkmer, rkmer, last_unknown);
            if (last_unknown >= K) {
                if (fkmer <= rkmer) {
                    // Is fwd
                    mers[fkmer]=0;
                } else {
                    // Is bwd
                    mers[rkmer]=0;
                }
            }
            ++s;
        }
    }
};
class KmerMapCounter : public  KMerFactory {
public:
    explicit KmerMapCounter(uint8_t k) : KMerFactory(k){}
    inline void count_all_kmers(const char * seq, std::unordered_map<uint64_t,uint32_t> &mers){
        // TODO: Adjust for when K is larger than what fits in uint64_t!
        last_unknown=0;
        fkmer=0;
        rkmer=0;
        auto s=seq;
        while (*s!='\0' and *s!='\n') {
            //fkmer: grows from the right (LSB)
            //rkmer: grows from the left (MSB)
            fillKBuf(*s, 0, fkmer, rkmer, last_unknown);
            if (last_unknown >= K) {
                if (fkmer <= rkmer) {
                    // Is fwd
                    auto m=mers.find(fkmer);
                    if (m!=mers.end()) ++m->second;
                } else {
                    // Is bwd
                    auto m=mers.find(rkmer);
                    if (m!=mers.end()) ++m->second;
                }
            }
            ++s;
        }
    }
};
class KmerVectorCreator : public  KMerFactory {
public:
    explicit KmerVectorCreator(uint8_t k) : KMerFactory(k){}
    inline std::vector<uint64_t> count_all_kmers(const char * seq){
        // TODO: Adjust for when K is larger than what fits in uint64_t!
        std::vector<uint64_t> v;
        last_unknown=0;
        fkmer=0;
        rkmer=0;
        auto s=seq;
        while (*s!='\0' and *s!='\n') {
            //fkmer: grows from the right (LSB)
            //rkmer: grows from the left (MSB)
            fillKBuf(*s, 0, fkmer, rkmer, last_unknown);
            if (last_unknown >= K) {
                if (fkmer <= rkmer) {
                    // Is fwd
                    v.emplace_back(fkmer);
                } else {
                    // Is bwd
                    v.emplace_back(rkmer);
                }
            }
            ++s;
        }
        return v;
    }
};

class UncoveredKmerCounter : public  KMerFactory {
public:
    explicit UncoveredKmerCounter(uint8_t k, const std::unordered_set<uint64_t> & _kset) : KMerFactory(k),kset(_kset){}
    inline uint64_t count_uncovered(const char * seq){
        // TODO: Adjust for when K is larger than what fits in uint64_t!

        last_unknown=0;
        fkmer=0;
        rkmer=0;
        auto s=seq;
        uint64_t uncovered=0;
        while (*s!='\0' and *s!='\n') {
            //fkmer: grows from the right (LSB)
            //rkmer: grows from the left (MSB)
            fillKBuf(*s, 0, fkmer, rkmer, last_unknown);
            if (last_unknown >= K) {
                if (fkmer <= rkmer) {
                    if (kset.count(fkmer)==0) ++uncovered;
                } else {
                    // Is bwd
                    if (kset.count(rkmer)==0) ++uncovered;
                }
            }
            ++s;
        }
        return uncovered;
    }
    const std::unordered_set<uint64_t> & kset;
};


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

size_t intersection_size_fast(const std::vector<bsg10xTag>& v1, const std::vector<bsg10xTag>& v2)
{
    size_t s=0;
    auto e1=v1.data()+v1.size();
    auto e2=v2.data()+v2.size();

    for (auto p1=v1.data(),p2=v2.data();p1<e1 and p2<e2;){
        if (*p1==*p2) {
            ++s;
            ++p1;
            ++p2;
        }
        else if (*p1<*p2) ++p1;
        else ++p2;
    }
    return s;
}

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

std::set<std::pair<sgNodeID_t, sgNodeID_t >> LinkageUntangler::get_HSPNPs(uint64_t min_size, float min_ci,
                                                                          float max_ci) {
    std::set<std::pair<sgNodeID_t, sgNodeID_t >> hspnps;
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
            else m = prev_fwl[0].dest;
            //Now evaluate coverage of the branches
            auto c1 = ws.kci.compute_compression_for_node(n, 1);
            if (std::isnan(c1) or c1<min_ci or c1>max_ci) continue;
            auto c2 = ws.kci.compute_compression_for_node(m, 1);
            if (std::isnan(c2) or c2<min_ci or c2>max_ci) continue;
#pragma omp critical(inserting_hspnps)
            {
                //hl<<(n<m ? n:m)<<" "<<(n<m ? m:n)<<std::endl;
                if (n < llabs(m)) hspnps.insert(std::make_pair(n, m));
                else hspnps.insert(std::make_pair(llabs(m), (m>0 ? n:-n)));
            }
        }
    }
    return hspnps;
}

void LinkageUntangler::select_nodes_by_HSPNPs(uint64_t min_size, float min_ci, float max_ci) {

    auto hspnps=get_HSPNPs(min_size,min_ci,max_ci);
    sglib::OutputLog() << "Selecting HSPNPs: " << hspnps.size() << " passed topology and CI" << std::endl;
    for (auto p:hspnps) {
        selected_nodes[llabs(p.first)] = true;
        selected_nodes[llabs(p.second)] = true;
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
    /*sglib::OutputLog()<<"filling orientation indexes"<<std::endl;
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
    sglib::OutputLog()<<"BW: "<<revc<<" ( "<<true_rev<<" - "<< false_rev<<" )"<<std::endl;*/
    std::map<std::pair<sgNodeID_t, sgNodeID_t>, uint64_t> lv;
    sglib::OutputLog()<<"collecting link votes across all paired libraries"<<std::endl;
    //use all libraries collect votes on each link
    auto rmi=0;
    for (auto &pm:ws.paired_read_mappers) {
        for (auto i = 1; i < pm.read_to_node.size(); i += 2) {
            sgNodeID_t n1 = pm.read_to_node[i];
            sgNodeID_t n2 = pm.read_to_node[i + 1];
            if (n1 == 0 or n2 == 0 or n1 == n2 or !selected_nodes[n1] or !selected_nodes[n2] ) continue;
            if (pm.read_direction_in_node[i]) n1=-n1;
            if (pm.read_direction_in_node[i+1]) n2=-n2;
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


LinkageDiGraph LinkageUntangler::make_tag_linkage(int min_reads, float end_perc) {

    //STEP 1 - identify candidates by simple tag-sharing.
    LinkageDiGraph ldg(ws.sg);

    //Step 1 - tag neighbours.

    sglib::OutputLog()<<"Getting tag neighbours"<<std::endl;
    auto pass_sharing=ws.linked_read_mappers[0].get_tag_neighbour_nodes(min_reads,selected_nodes);

    sglib::OutputLog()<<"Node pairs with more than "<<min_reads<<" shared tags: "<<pass_sharing.size()<<std::endl;

    //STEP 2 - confirm directionality

    //2.a create link direction counts:
    std::map<std::pair<sgNodeID_t, sgNodeID_t>, uint64_t> lv;
    sglib::OutputLog()<<"collecting link votes across all paired libraries"<<std::endl;
    //use all libraries collect votes on each link
    auto rmi=0;
    for (auto &pm:ws.paired_read_mappers) {
        for (auto i = 1; i < pm.read_to_node.size(); i += 2) {
            sgNodeID_t n1 = pm.read_to_node[i];
            sgNodeID_t n2 = pm.read_to_node[i + 1];
            if (n1 == 0 or n2 == 0 or n1 == n2 or !selected_nodes[n1] or !selected_nodes[n2] ) continue;
            if (pm.read_direction_in_node[i]) n1=-n1;
            if (pm.read_direction_in_node[i+1]) n2=-n2;
            if (llabs(n1) > llabs(n2)) std::swap(n1,n2);
            ++lv[std::make_pair(n1, n2)];
        }
        ++rmi;
    }
    std::set<std::pair<sgNodeID_t ,sgNodeID_t >> used;
    for (auto p:pass_sharing) {
        auto bf=lv[std::make_pair(-p.first,p.second)];
        auto bb=lv[std::make_pair(-p.first,-p.second)];
        auto ff=lv[std::make_pair(p.first,p.second)];
        auto fb=lv[std::make_pair(p.first,-p.second)];
        auto total=bf+bb+ff+fb;
        float bfp=((float) bf)/total;
        float bbp=((float) bb)/total;
        float ffp=((float) ff)/total;
        float fbp=((float) fb)/total;
        if (bf>=3 and bfp>=.75) {
            ldg.add_link(-p.first,p.second,0);
            used.insert(p);
        }
        else if (bb>=3 and bbp>=.75) {
            ldg.add_link(-p.first,-p.second,0);
            used.insert(p);
        }
        else if (ff>=3 and ffp>=.75) {
            ldg.add_link(p.first,p.second,0);
            used.insert(p);
        }
        else if (fb>=3 and fbp>=.75) {
            ldg.add_link(p.first,-p.second,0);
            used.insert(p);
        }
        /*std::cout<<"Evaluating connection between "<<p.first<<" and "<<p.second<<": "
                <<lv[std::make_pair(-p.first,p.second)]<<" "
                <<lv[std::make_pair(-p.first,-p.second)]<<" "
                <<lv[std::make_pair(p.first,p.second)]<<" "
                <<lv[std::make_pair(p.first,-p.second)]<<std::endl;*/
    }

    //STEP 3 - Looking at disconnected ends on 1-0 and N-0 nodes
    std::vector<sgNodeID_t> one_end_only;
    uint64_t disc=0,ldisc=0,single=0,lsingle=0,both=0,lboth=0;
    for (sgNodeID_t n=1;n<ws.sg.nodes.size();++n) {
        if (!selected_nodes[n]) continue;
        auto blc=ldg.get_bw_links(n).size();
        auto flc=ldg.get_fw_links(n).size();

        if (blc==0 and flc==0){
            ++disc;
            if (ws.sg.nodes[n].sequence.size()>2000) ++ldisc;
        }
        else if (blc==0 or flc==0){
            if (blc==0) one_end_only.push_back(-n);
            else one_end_only.push_back(n);
            ++single;
            if (ws.sg.nodes[n].sequence.size()>2000) ++lsingle;
        } else {
            ++both;
            if (ws.sg.nodes[n].sequence.size()>2000) ++lboth;
        }
    }
    /*sglib::OutputLog()<<both<<" nodes with both-sides linkage ( "<<lboth<<" >2kbp )"<<std::endl;
    sglib::OutputLog()<<single<<" nodes with one-side linkage ( "<<lsingle<<" >2kbp )"<<std::endl;
    sglib::OutputLog()<<disc<<" nodes without linkage ( "<<ldisc<<" >2kbp )"<<std::endl;*/
    ldg.report_connectivity();
    sglib::OutputLog()<<"Attempting single-side reconnection through topology"<<std::endl;
    auto tldg=make_topology_linkage(30);
#pragma omp parallel for
    for (auto i=0; i<one_end_only.size();++i){
        auto n=one_end_only[i];
        //first look for the topology connection.
        for (auto tfnl:tldg.get_fw_links(n)){
            std::pair<sgNodeID_t, sgNodeID_t> pair;
            pair.first=llabs(n);
            pair.second=llabs(tfnl.dest);
            if (pair.first>pair.second) std::swap(pair.first,pair.second);
            for (auto ps:pass_sharing) if (ps==pair) {
#pragma omp critical (add_topo_link)
                ldg.add_link(tfnl.source,tfnl.dest,0);
            }
        }
    }
    ldg.report_connectivity();

    /*sglib::OutputLog()<<"Evaluating tag imbalance"<<std::endl;
    for (auto p:pass_sharing) {

        auto n1 = p.first;
        auto n2 = p.second;
        std::set<bsg10xTag> shared_tags;
        std::set_intersection(node_tags[n1].begin(), node_tags[n1].end(), node_tags[n2].begin(), node_tags[n2].end(),
                              std::inserter(shared_tags, shared_tags.end()));
        uint64_t n1_front_in = 0, n1_front_total = 0, n1_back_in = 0, n1_back_total = 0;
        uint64_t n2_front_in = 0, n2_front_total = 0, n2_back_in = 0, n2_back_total = 0;
        uint64_t n1first30point = ws.sg.nodes[n1].sequence.size() * end_perc;
        uint64_t n1last30point = ws.sg.nodes[n1].sequence.size() * (1 - end_perc);
        std::set<bsg10xTag> t1f,t1b,t2f,t2b,t1ft,t1bt,t2ft,t2bt;
        for (auto rm:ws.linked_read_mappers[0].reads_in_node[n1]) {
            if (rm.first_pos < n1first30point) {
                ++n1_front_total;
                t1ft.insert(ws.linked_read_datastores[0].get_read_tag(rm.read_id));
                if (shared_tags.count(ws.linked_read_datastores[0].get_read_tag(rm.read_id)) > 0) {
                    ++n1_front_in;
                    t1f.insert(ws.linked_read_datastores[0].get_read_tag(rm.read_id));
                }
            }
            if (rm.last_pos > n1last30point) {
                ++n1_back_total;
                t1bt.insert(ws.linked_read_datastores[0].get_read_tag(rm.read_id));
                if (shared_tags.count(ws.linked_read_datastores[0].get_read_tag(rm.read_id)) > 0) {
                    ++n1_back_in;
                    t1b.insert(ws.linked_read_datastores[0].get_read_tag(rm.read_id));
                }
            }
        }
        auto n1f = (100.0 * n1_front_in / n1_front_total);
        auto n1b = (100.0 * n1_back_in / n1_back_total);
        uint64_t n2first30point = ws.sg.nodes[n2].sequence.size() * end_perc;
        uint64_t n2last30point = ws.sg.nodes[n2].sequence.size() * (1 - end_perc);
        for (auto rm:ws.linked_read_mappers[0].reads_in_node[n2]) {
            if (rm.first_pos < n2first30point) {
                ++n2_front_total;
                t2ft.insert(ws.linked_read_datastores[0].get_read_tag(rm.read_id));
                if (shared_tags.count(ws.linked_read_datastores[0].get_read_tag(rm.read_id)) > 0) {
                    ++n2_front_in;
                    t2f.insert(ws.linked_read_datastores[0].get_read_tag(rm.read_id));
                }
            }
            if (rm.last_pos > n2last30point) {
                ++n2_back_total;
                t2bt.insert(ws.linked_read_datastores[0].get_read_tag(rm.read_id));
                if (shared_tags.count(ws.linked_read_datastores[0].get_read_tag(rm.read_id)) > 0) {
                    ++n2_back_in;
                    t2b.insert(ws.linked_read_datastores[0].get_read_tag(rm.read_id));
                }
            }
        }
        auto n2f = (100.0 * n2_front_in / n2_front_total);
        auto n2b = (100.0 * n2_back_in / n2_back_total);
        if ( (ws.sg.nodes[llabs(n1)].sequence.size()>10000 and ws.sg.nodes[llabs(n2)].sequence.size()>10000) ){
            std::cout<<"connection between "<<n1<<" and "<<n2<<" with "<<shared_tags.size()<<" tags: "<<n1f<<"("<<t1f.size()<<"):"<< n1b <<"("<<t1b.size()<<") <-> "<<n2f<<"("<<t2f.size()<<"):"<< n2b <<"("<<t2b.size()<<")"<<std::endl;
            std::cout<<"F<->F: "<<intersection_size(t1f,t2f)<<" / "<<t1ft.size()<<":"<<t2ft.size();
            std::cout<<"  F<->B: "<<intersection_size(t1f,t2b)<<" / "<<t1ft.size()<<":"<<t2bt.size();
            std::cout<<"  B<->F: "<<intersection_size(t1b,t2f)<<" / "<<t1bt.size()<<":"<<t2ft.size();
            std::cout<<"  B<->B: "<<intersection_size(t1b,t2b)<<" / "<<t1bt.size()<<":"<<t2bt.size()<<std::endl;
        }
        if (fabs(2 * (n1f - n1b) / (n1f + n1b)) > .1 and fabs(2 * (n2f - n2b) / (n2f + n2b)) > .1) {
#pragma omp critical
            ++linked;
            ldg.add_link((n1f > n1b ? n1 : -n1), (n2f > n2b ? n2 : -n2), 0);
        }
    }
    sglib::OutputLog()<<"Links created (passing tag imbalance): "<<linked<<std::endl;*/
    return ldg;
}

LinkageDiGraph LinkageUntangler::filter_linkage_to_hspnp_duos(uint64_t min_size, float min_ci, float max_ci,
                                                              const LinkageDiGraph &ldg_old) {
    std::unordered_map<sgNodeID_t,sgNodeID_t> node_to_parallel;
    //1- get all hspnps -> create a map of parallels
    LinkageDiGraph ldg_new(ws.sg);
    auto hspnps=get_HSPNPs(min_size,min_ci,max_ci);
    for (auto h:hspnps) {
        node_to_parallel[h.first]=h.second;
        node_to_parallel[-h.first]=-h.second;
        node_to_parallel[h.second]=h.first;
        node_to_parallel[-h.second]=-h.first;
    }
    //2- hspnp -> look for links in one direction from one of the nodes, and same direction for the other
    for (auto h:hspnps){
        auto hr=h;
        hr.first=-hr.first;
        hr.second=-hr.second;
        for (auto hspnp:{h,hr}) {
            auto n1fs = ldg_old.get_fw_links(hspnp.first);
            auto n2fs = ldg_old.get_fw_links(hspnp.second);
            for (auto n1f:n1fs) {
                for (auto n2f:n2fs) {
                    if (node_to_parallel.count(n1f.dest) and node_to_parallel[n1f.dest] == n2f.dest) {
                        // if links are to parts of the same node -> introduce linkage on newldg.
                        ldg_new.add_link(-hspnp.first, n1f.dest, 0);
                        ldg_new.add_link(-hspnp.second, n2f.dest, 0);
                    }
                }
            }
        }
    }
    return ldg_new;

}

void LinkageUntangler::expand_trivial_repeats(const LinkageDiGraph & ldg) {
    uint64_t aa=0,ab=0;
    for (auto n=1;n<ws.sg.nodes.size();++n) {
        if (ws.sg.nodes[n].status == sgNodeDeleted) continue;
        //check node is 2-2
        auto bwl=ws.sg.get_bw_links(n);
        if (bwl.size()!=2) continue;
        auto fwl=ws.sg.get_fw_links(n);
        if (fwl.size()!=2) continue;
        auto p1=-bwl[0].dest;
        auto p2=-bwl[1].dest;
        auto n1=fwl[0].dest;
        auto n2=fwl[1].dest;
        //check bw nodes have only one fw, is one of the fws and not the same
        auto p1ll=ldg.get_fw_links(p1);
        if (p1ll.size()!=1) continue;
        auto p2ll=ldg.get_fw_links(p2);
        if (p2ll.size()!=1) continue;
        if (p1ll[0].dest==n1 and p2ll[0].dest==n2){
            ws.sg.expand_node(n,{{p1},{p2}},{{n1},{n2}});
            ++aa;
        }
        else if (p2ll[0].dest==n1 and p1ll[0].dest==n2) {
            ws.sg.expand_node(n,{{p1},{p2}},{{n2},{n1}});
            ++ab;
        }
        else continue;
    }
    sglib::OutputLog()<<"Repeat expansion: AA:"<<aa<<"  AB:"<<ab<<std::endl;
}

void LinkageUntangler::expand_linear_regions(const LinkageDiGraph & ldg) {
    sglib::OutputLog()<<"Starting linear region expansion..."<<std::endl;
    //sglib::OutputLog()<<"Looking for \"lines\"..."<<std::endl;
    auto lines=ldg.get_all_lines(2);
    sglib::OutputLog()<<"Creating tag sets for "<<lines.size()<<" linear regions"<<std::endl;
    //sglib::OutputLog()<<"TODO: now use tags and LMPs to find paths between elements in the line"<<std::endl;
    //sglib::OutputLog()<<"USING ONLY 10 lines as a test"<<std::endl;
    //lines.resize(10);
    //---------------------------------Step 1: get tagsets for lines.
    std::vector<std::set<bsg10xTag>> linetagsets;
    linetagsets.reserve(lines.size());
    BufferedTagKmerizer btk(ws.linked_read_datastores[0],31,100000,1000);
    for (auto l:lines){
        //sglib::OutputLog()<<"Analising line: ";
        //for (auto &ln:l) std::cout<<"seq"<<llabs(ln)<<", ";
        //for (auto &ln:l) std::cout<<ln<<" ";
        //std::cout<<std::endl;
        std::map<bsg10xTag ,std::pair<uint32_t , uint32_t >> tagcounts; //tag -> nodes, reads
        for (auto &ln:l) {
            std::map<bsg10xTag ,uint32_t> ntagcounts;
            for (auto rm:ws.linked_read_mappers[0].reads_in_node[llabs(ln)]){
                auto tag=ws.linked_read_datastores[0].get_read_tag(rm.read_id);
                ++ntagcounts[tag];
            }
            for (auto ntc:ntagcounts) {
                ++tagcounts[ntc.first].first;
                tagcounts[ntc.first].second+=ntc.second;
            }
        }
        std::map<bsg10xTag ,std::pair<uint32_t , uint32_t >> tagtotals;
        std::set<bsg10xTag> lineTagSet;
        for (auto tc:tagcounts) {
            auto tag=tc.first;
            auto reads=ws.linked_read_datastores[0].get_tag_reads(tc.first);
            std::set<sgNodeID_t> nodes;
            for (auto r:reads) nodes.insert(ws.linked_read_mappers[0].read_to_node[r]);
            tagtotals[tag].first=nodes.size()-nodes.count(0);
            tagtotals[tag].second=reads.size();
            if (tc.second.first>1 and reads.size()<3000) lineTagSet.insert(tc.first);
        }
        linetagsets.push_back(lineTagSet);
        if (linetagsets.size()%100==0) std::cout<<"."<<std::flush;
    }
    std::cout<<std::endl;
    sglib::OutputLog()<<"Creating path collections to be evaluated for "<<lines.size()<<" linear regions"<<std::endl;
    std::vector<std::vector<std::vector<SequenceGraphPath>>> alternatives;
    uint64_t total_paths=0,found=0,evaluated=0;
    alternatives.reserve(lines.size());
    for (auto l:lines) {
        alternatives.emplace_back();
        for (auto i = 0; i < l.size() - 1; ++i) {
            evaluated++;
            auto from = l[i];
            auto to = l[i + 1];
            auto paths = ws.sg.find_all_paths_between(from, to, 400000, 20);
            if (paths.size()>0) found++;
            alternatives.back().emplace_back(paths);
            total_paths+=paths.size();
            //sglib::OutputLog() << paths.size() << " paths to go from " << from << " to " << to << std::endl;
        }
        if (alternatives.size()%100==0) std::cout<<"."<<std::flush;
    }
    std::cout<<std::endl;
    sglib::OutputLog()<<"Junctions with possible paths: "<<found<<" / "<<evaluated<<std::endl;
    sglib::OutputLog()<<"Total paths to evaluate: "<<total_paths<<std::endl;
    //Now use a function that only counts coverage on the set of kmers from all paths collections for each line
    //kmer_coverage_in_tagreads(&std::map<kmer, coverage> (init at 0), std::set<tag> linetagset);

    std::cout << "creating and populating the maps as of now" << std::endl;
    std::vector<std::unordered_map<uint64_t, uint32_t>> linekmercoverages;
    linekmercoverages.resize(lines.size());
#pragma omp parallel
    {
        KmerMapCounter km_count(31);
        KmerMapCreator km_create(31);
        BufferedLRSequenceGetter blrsg(ws.linked_read_datastores[0], 200000, 1000);
        std::unordered_map<uint64_t, uint32_t> kmercoverages;
        uint64_t done=0;
#pragma omp for schedule(static, 100)
        for (auto i = 0; i < lines.size(); ++i) {
            //map with all kmers of paths to be evaluated
            kmercoverages.clear();
//            size_t t=0;
//            for (auto &alts:alternatives[i]) {
//                for (auto &a:alts) {
//                    for (auto n:a.nodes) t += ws.sg.nodes[llabs(n)].sequence.size();
//                }
//            }
//            kmercoverages.reserve(t);
            for (auto &alts:alternatives[i]) {
                for (auto &a:alts) {
                    for (auto n:a.nodes) {
                        km_create.create_all_kmers(ws.sg.nodes[llabs(n)].sequence.c_str(), kmercoverages);
                    }
                }
            }
            for (auto &t:linetagsets[i]) {
                for (auto rid:ws.linked_read_datastores[0].get_tag_reads(t)) {
                    km_count.count_all_kmers(blrsg.get_read_sequence(rid), kmercoverages);
                }
            }
#pragma omp critical
            linekmercoverages[i]=kmercoverages;
            ++done;
            if (done % 100 == 0) std::cout << "." << std::flush;
            //count from the tag's reads
            //btk.get_tag_kmers()
        }
    }
    std::cout<<"DONE"<<std::endl;
    sglib::OutputLog()<<"evaluating alternative paths between each pair of adjacent nodes"<<std::endl;
    KmerVectorCreator kvc(31);
    uint64_t solved=0,none_covered=0,too_many_covered=0,no_paths=0;
    GraphEditor ged(ws);
    std::vector<SequenceGraphPath> sols;
    for (auto i=0;i<lines.size();++i){
        for (auto ia=0;ia<alternatives[i].size();++ia){
            int best=-1;
            bool too_many=false;
            for (auto j=0;j<alternatives[i][ia].size();++j){
                uint64_t missed=0;
                for (auto n:alternatives[i][ia][j].nodes) {
                    for (auto x:kvc.count_all_kmers(ws.sg.nodes[llabs(n)].sequence.c_str())) {
                        if (linekmercoverages[i][x] < 8) ++missed;//TODO: maybe ask for more than 1 read coverage?
                    }
                }
                if (missed==0){
                    if (best==-1) {
                        best = j;
                    }
                    else {
                        too_many=true;
                        best = -1;
                        break;
                    }
                }
            }
            //std::cout<<"Solution for line "<<i<<" jump #"<<ia<<": "<<best<<" / "<<alternatives[i][ia].size() <<std::endl;
            if (best!=-1){
                for (auto n:alternatives[i][ia][best].nodes) {
                    if (selected_nodes[llabs(n)]){
                        best=-1;
                        break;
                    }
                }
            }
            if (best==-1) {
                if (alternatives[i][ia].empty()) ++no_paths;
                else if (too_many) ++too_many_covered;
                else ++none_covered;
            }
            else {
                ++solved;
                sols.emplace_back(ws.sg);
                sols.back().nodes.emplace_back(lines[i][ia]);
                for (auto n:alternatives[i][ia][best].nodes) sols.back().nodes.emplace_back(n);
                sols.back().nodes.emplace_back(lines[i][ia+1]);
            }
        }
    }
    std::cout<<"Solved: "<<solved<<"  Too many covered paths: "<<too_many_covered<<"  No covered paths: "<<none_covered<<"  No paths found: "<<no_paths<<std::endl;
    sglib::OutputLog()<<"Applying solutions in the graph"<<std::endl;
    uint64_t applied=0;
    for (auto s:sols) {
        if (ged.detach_path(s)) ++applied;
    }
    sglib::OutputLog()<<applied<<" solutions applied"<<std::endl;
}

void LinkageUntangler::linear_regions_tag_local_assembly(const LinkageDiGraph & ldg, uint8_t k, int min_cvg, int max_lines, uint64_t min_nodes, uint64_t min_total_size){
    sglib::OutputLog()<<"Starting linear region tag local assemblies..."<<std::endl;
    auto lines=ldg.get_all_lines(min_nodes, min_total_size);
    if (max_lines>0) {
        sglib::OutputLog()<<"USING ONLY "<<max_lines<< " lines as a test"<<std::endl;
        lines.resize(max_lines);
        //to filter to specific node, keep only the line that has that node:
        //auto old_lines=lines;
        //lines.clear();
        //for (auto l:lines) for (auto n:l) if (llabs(n)==TARGET_NODE_ID) lines.push_back(l);
    }

    sglib::OutputLog()<<"Creating tag sets for "<<lines.size()<<" linear regions"<<std::endl;
    //---------------------------------Step 1: get tagsets for lines.
    std::vector<std::set<bsg10xTag>> linetagsets;
    linetagsets.reserve(lines.size());
    for (auto l:lines){
        std::map<bsg10xTag ,std::pair<uint32_t , uint32_t >> tagcounts; //tag -> nodes, reads
        for (auto &ln:l) {
            std::map<bsg10xTag ,uint32_t> ntagcounts;
            for (auto rm:ws.linked_read_mappers[0].reads_in_node[llabs(ln)]){
                auto tag=ws.linked_read_datastores[0].get_read_tag(rm.read_id);
                ++ntagcounts[tag];
            }
            for (auto ntc:ntagcounts) {
                ++tagcounts[ntc.first].first;
                tagcounts[ntc.first].second+=ntc.second;
            }
        }
        std::map<bsg10xTag ,std::pair<uint32_t , uint32_t >> tagtotals;
        std::set<bsg10xTag> lineTagSet;
        for (auto tc:tagcounts) {
            auto tag=tc.first;
            auto reads=ws.linked_read_datastores[0].get_tag_reads(tc.first);
            std::set<sgNodeID_t> nodes;
            for (auto r:reads) nodes.insert(ws.linked_read_mappers[0].read_to_node[r]);
            tagtotals[tag].first=nodes.size()-nodes.count(0);
            tagtotals[tag].second=reads.size();
            if (tc.second.first>1 and reads.size()<3000) lineTagSet.insert(tc.first);
        }
        linetagsets.push_back(lineTagSet);
        if (linetagsets.size()%100==0) std::cout<<"."<<std::flush;
    }
    std::cout<<std::endl;
    uint64_t jc=0;
    for (auto &l:lines) jc+=l.size()-1;

    std::vector<SequenceGraphPath> sols;
#pragma omp parallel
    {
        BufferedLRSequenceGetter blrsg(ws.linked_read_datastores[0], 200000, 1000);
        std::vector<SequenceGraphPath> tsols;
        uint64_t donelines=0;
#pragma omp for schedule(dynamic,1)
        for (auto i=0; i<lines.size(); ++i){
            auto ltkmers128 = ws.linked_read_datastores[0].get_tags_kmers128(k, min_cvg, linetagsets[i], blrsg);
            std::cout << "creating DBG for line #" << i << std::endl;
            SequenceGraph dbg;
            GraphMaker gm(dbg);
            gm.new_graph_from_kmerset_trivial128(ltkmers128, k);
            //dbg.write_to_gfa("local_dbg_" + std::to_string(i) + ".gfa");
            //gruesome tip clipping:
            std::cout << "Starting gruesome tip clipping" << std::endl;
            std::set<sgNodeID_t> to_delete;
            for (sgNodeID_t n = 1; n < dbg.nodes.size(); ++n) {
                if (dbg.nodes[n].status == sgNodeDeleted) continue;
                if (dbg.nodes[n].sequence.size() > 200) continue;
                //std::cout<<"Evaluating seq"<<n<<": ";
                auto fwl = dbg.get_fw_links(n);
                auto bwl = dbg.get_bw_links(n);
                //std::cout<<" fwl: "<<fwl.size()<<"  bwl: "<<bwl.size();
                if (fwl.size() == 1 and bwl.size() == 0) {
                    //std::cout<<"  bwl for "<<fwl[0].dest<<": "<<dbg.get_bw_links(fwl[0].dest).size();
                    if (dbg.get_bw_links(fwl[0].dest).size() == 2) {
                        to_delete.insert(n);
                        //std::cout<<" D"<<std::endl;
                    }
                }
                if (fwl.size() == 0 and bwl.size() == 1) {
                    //std::cout<<"  fwl for "<<-bwl[0].dest<<": "<<dbg.get_fw_links(-bwl[0].dest).size();
                    if (dbg.get_fw_links(-bwl[0].dest).size() == 2) {
                        to_delete.insert(n);
                        //std::cout<<" D"<<std::endl;
                    }
                }
                if (fwl.size() == 0 and bwl.size()==0) to_delete.insert(n);
                //std::cout<<std::endl;
            }
            std::cout << "Nodes to delete: " << to_delete.size() << std::endl;
            for (auto n:to_delete) dbg.remove_node(n);
            auto utc = dbg.join_all_unitigs();
            std::cout << "Joined unitigs: " << utc << std::endl;
            dbg.write_to_gfa("local_dbg_" + std::to_string(i) + "_tipclipped.gfa");
            std::ofstream anchf("local_dbg_" + std::to_string(i) + "_anchors.fasta");
            for (auto n:lines[i]){
                anchf<<">seq"<<llabs(n)<<std::endl;
                anchf<<ws.sg.nodes[llabs(n)].sequence<<std::endl;

            }
            if (++donelines%100==0) std::cout<<"."<<std::flush;
        }
    }
}

void LinkageUntangler::expand_linear_regions_skating(const LinkageDiGraph & ldg, int max_lines) {
    sglib::OutputLog()<<"Starting linear region consolidation via skating with line tag collection..."<<std::endl;
    auto lines=ldg.get_all_lines(2);
    if (max_lines>0) {
        sglib::OutputLog()<<"USING ONLY "<<max_lines<< " lines as a test"<<std::endl;
        lines.resize(max_lines);
    }

    sglib::OutputLog()<<"Creating tag sets for "<<lines.size()<<" linear regions"<<std::endl;
    //---------------------------------Step 1: get tagsets for lines.
    std::vector<std::set<bsg10xTag>> linetagsets;
    linetagsets.reserve(lines.size());
    BufferedTagKmerizer btk(ws.linked_read_datastores[0],31,100000,1000);
    for (auto l:lines){
        //sglib::OutputLog()<<"Analising line: ";
        //for (auto &ln:l) std::cout<<"seq"<<llabs(ln)<<", ";
        //for (auto &ln:l) std::cout<<ln<<" ";
        //std::cout<<std::endl;
        std::map<bsg10xTag ,std::pair<uint32_t , uint32_t >> tagcounts; //tag -> nodes, reads
        for (auto &ln:l) {
            std::map<bsg10xTag ,uint32_t> ntagcounts;
            for (auto rm:ws.linked_read_mappers[0].reads_in_node[llabs(ln)]){
                auto tag=ws.linked_read_datastores[0].get_read_tag(rm.read_id);
                ++ntagcounts[tag];
            }
            for (auto ntc:ntagcounts) {
                ++tagcounts[ntc.first].first;
                tagcounts[ntc.first].second+=ntc.second;
            }
        }
        std::map<bsg10xTag ,std::pair<uint32_t , uint32_t >> tagtotals;
        std::set<bsg10xTag> lineTagSet;
        for (auto tc:tagcounts) {
            auto tag=tc.first;
            auto reads=ws.linked_read_datastores[0].get_tag_reads(tc.first);
            std::set<sgNodeID_t> nodes;
            for (auto r:reads) nodes.insert(ws.linked_read_mappers[0].read_to_node[r]);
            tagtotals[tag].first=nodes.size()-nodes.count(0);
            tagtotals[tag].second=reads.size();
            if (tc.second.first>1 and reads.size()<3000) lineTagSet.insert(tc.first);
        }
        linetagsets.push_back(lineTagSet);
        if (linetagsets.size()%100==0) std::cout<<"."<<std::flush;
    }
    std::cout<<std::endl;
    uint64_t jc=0;
    for (auto &l:lines) jc+=l.size()-1;
    sglib::OutputLog()<<"Skating across "<<jc<<" junctions in "<<lines.size()<<" linear regions"<<std::endl;

    std::vector<SequenceGraphPath> sols;
#pragma omp parallel
    {
        BufferedLRSequenceGetter blrsg(ws.linked_read_datastores[0], 200000, 1000);
        std::vector<SequenceGraphPath> tsols;
        uint64_t donelines=0;
#pragma omp for schedule(dynamic,1)
        for (auto i=0; i<lines.size(); ++i){
            //std::cout<<"Creating kmer set for line"<<i<<" from tags"<<std::endl;
            auto ltkmers=ws.linked_read_datastores[0].get_tags_kmers(31,3,linetagsets[i],blrsg);
            UncoveredKmerCounter ukc(31,ltkmers);
            //std::cout<<"Evaluating paths for "<<lines[i].size()-1<<" junctions"<<std::endl;
            for (auto j=0;j<lines[i].size()-1;++j){
                auto from=lines[i][j];
                auto to=lines[i][j+1];
                //std::cout<<std::endl<<std::endl<<"Junction #"<<j+1<<" from "<<from<<" to "<<to<<std::endl;
                std::vector<std::vector<sgNodeID_t>> skated_paths;
                skated_paths.push_back({from});
                int max_nodes=50;
                while (--max_nodes and not skated_paths.empty()){
                    //std::cout<<std::endl<<"expansion round starting with "<<skated_paths.size()<<" paths "<<std::endl;
                    auto old_skated=skated_paths;
                    skated_paths.clear();
                    bool loop=false,crosstalk=false;
                    for (auto p:old_skated) {
                        if (p.back()==to) {
                            skated_paths.push_back(p);
                            continue;
                        }
                        //std::cout<<" expanding fw from node "<<p.back()<<std::endl;
                        for (auto fwl:ws.sg.get_fw_links(p.back())) {
                            //std::cout<<"  considering fwl to "<<fwl.dest<<std::endl;
                            if (std::count(p.begin(),p.end(),fwl.dest)>0 or std::count(p.begin(),p.end(),-fwl.dest)>0){
                                loop=true;
                                //std::cout<<"loop detected, aborting junction analysis"<<std::endl;
                                break;
                            }

                            auto u=ukc.count_uncovered(ws.sg.nodes[llabs(fwl.dest)].sequence.c_str());
                            //std::cout<<"  Uncovered kmers in "<<fwl.dest<<" ("<<ws.sg.nodes[llabs(fwl.dest)].sequence.size()<<" bp): "
                            //                                                                                                <<u<<std::endl;
                            if ( u == 0) {
                                //check for a path that reaches a selected node that is not connected here
                                if (selected_nodes[llabs(fwl.dest)] and fwl.dest!=to) {
                                    crosstalk=true;
                                    break;
                                }
                                //std::cout<<"  path can continue in node"<<fwl.dest<<std::endl;
                                skated_paths.push_back(p);
                                skated_paths.back().push_back(fwl.dest);
                            }
                        }
                    }
                    if (loop or crosstalk) {
                        skated_paths.clear();
                        break;
                    }
                }
                uint64_t complete=0,incomplete=0;
                for (auto p:skated_paths) {
                    if (p.back()==to) ++complete;
                    else ++incomplete;
                }
                if (complete==1 and incomplete==0) tsols.emplace_back(SequenceGraphPath(ws.sg,skated_paths[0]));
                //std::cout<<"Skating line #"<<i+1<<" junction #"<<j+1<<" produced "<<complete<<" complete paths and "<<incomplete<<" possibly incomplete paths"<<std::endl;
            }
            if (++donelines%100==0) std::cout<<"."<<std::flush;
        }
#pragma omp critical
        sols.insert(sols.end(),tsols.begin(),tsols.end());
    }
    sglib::OutputLog()<<"Applying "<<sols.size()<<" solutions in the graph"<<std::endl;
    GraphEditor ged(ws);
    uint64_t applied=0;
    for (auto s:sols) {
        if (ged.detach_path(s)) ++applied;
    }
    sglib::OutputLog()<<applied<<" solutions applied"<<std::endl;
}