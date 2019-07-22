//
// Created by Bernardo Clavijo (EI) on 28/05/2018.
//

#include "LinkageUntangler.hpp"
#include "GraphEditor.hpp"
#include "GraphMaker.hpp"
#include "LinkageMaker.hpp"
#include <atomic>

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
            fillKBuf(*s, fkmer, rkmer, last_unknown);
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
            fillKBuf(*s, fkmer, rkmer, last_unknown);
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
            fillKBuf(*s, fkmer, rkmer, last_unknown);
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
            fillKBuf(*s, fkmer, rkmer, last_unknown);
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

//======= NODE SELECTION METHODS =======

void LinkageUntangler::deselect_all() {
    selected_nodes.clear();
    selected_nodes.resize(dg.sdg.nodes.size(),false);
}

void LinkageUntangler::select_all() {
    selected_nodes.clear();
    selected_nodes.resize(dg.sdg.nodes.size(),true);
}

void LinkageUntangler::report_selection() {
    uint64_t total_bp=0,total_count=0,selected_bp=0,selected_count=0;
    for (auto n=1;n<dg.sdg.nodes.size();++n) {
        if (dg.sdg.nodes[n].status == NodeStatus::Deleted) continue;
        total_bp+=dg.sdg.nodes[n].sequence.size();
        ++total_count;
        if (selected_nodes[n]) {
            selected_bp += dg.sdg.nodes[n].sequence.size();
            ++selected_count;
        }
    }
        sdglib::OutputLog()<< "Current selection: "<<selected_count<<" / "<<total_count<<" nodes  with  "<<selected_bp<<" / "<<total_bp<<" bp"<<std::endl;

}
void LinkageUntangler::select_by_size( uint64_t min_size, uint64_t max_size) {
    {
        for (auto n=1;n<dg.sdg.nodes.size();++n) {
            if (dg.sdg.nodes[n].status==NodeStatus::Deleted) continue;
            if (dg.sdg.nodes[n].sequence.size() >= min_size and
                (max_size==0 or dg.sdg.nodes[n].sequence.size() <= max_size))
                selected_nodes[n]=true;
        }
    }
}

//void LinkageUntangler::select_nodes_by_size_and_ci( uint64_t min_size, float min_ci, float max_ci) {
//    std::vector<sgNodeID_t> nodes;
//    sdglib::OutputLog()<<"LU selecting nodes by size and ci: size >= " << min_size << " bp  |  " << min_ci << "<= CI <=" << max_ci <<std::endl;
//    uint64_t deleted=0,small=0,cifail=0,selected=0;
//#pragma omp parallel
//    {
//#pragma omp for schedule(static, 100) reduction(+:deleted,small,cifail,selected)
//        for (auto n=1;n<ws.sdg.nodes.size();++n) {
//            if (ws.sdg.nodes[n].status==NodeStatus::Deleted) { ++deleted; continue; }
//            if (ws.sdg.nodes[n].sequence.size() < min_size) { ++small; continue; }
//            auto ci = ws.kci.compute_compression_for_node(n, 1);
//            if (std::isnan(ci) or ci < min_ci or ci > max_ci) { ++cifail; continue;}
//            #pragma omp critical(collect_selected_nodes)
//            selected_nodes[n]=true;
//            ++selected;
//        }
//    }
//    sdglib::OutputLog()<<deleted<<" deleted, "<<small<<" small, "<<cifail<<" wrong CI and "<<selected<<" selected nodes."<<std::endl;
//
//}

//std::set<std::pair<sgNodeID_t, sgNodeID_t >> LinkageUntangler::get_HSPNPs(uint64_t min_size, float min_ci,
//                                                                          float max_ci) {
//    std::set<std::pair<sgNodeID_t, sgNodeID_t >> hspnps;
//#pragma omp parallel for schedule(static, 100)
//    for (sgNodeID_t n = 1; n < ws.sdg.nodes.size(); ++n) {
//        if (ws.sdg.nodes[n].status == NodeStatus::Deleted) continue;
//        if (ws.sdg.nodes[n].sequence.size() < min_size) continue;
//        //FW check
//        auto fwl = ws.sdg.get_fw_links(n);
//        if (fwl.size() != 1) continue;
//        auto post = fwl[0].dest;
//        auto post_bwl = ws.sdg.get_bw_links(post);
//        if (post_bwl.size() != 2) continue;
//        if (llabs(post_bwl[0].dest)==llabs(post_bwl[1].dest))continue;
//        //BW check
//        auto bwl = ws.sdg.get_bw_links(n);
//        if (bwl.size() != 1) continue;
//        auto prev = bwl[0].dest;
//        auto prev_fwl = ws.sdg.get_bw_links(prev);
//        if (prev_fwl.size() != 2) continue;
//
//        if ((prev_fwl[0].dest == -post_bwl[0].dest and prev_fwl[1].dest == -post_bwl[1].dest)
//            or (prev_fwl[1].dest == -post_bwl[0].dest and prev_fwl[0].dest == -post_bwl[1].dest)) {
//            sgNodeID_t m;
//            if (llabs(prev_fwl[0].dest) != n and llabs(prev_fwl[1].dest) != n) std::cout<<"Error! cant find N in prev!"<<std::endl;
//            if (llabs(prev_fwl[0].dest) == n) m = llabs(prev_fwl[1].dest);
//            else m = prev_fwl[0].dest;
//            //Now evaluate coverage of the branches
//            auto c1 = ws.kci.compute_compression_for_node(n, 1);
//            if (std::isnan(c1) or c1<min_ci or c1>max_ci) continue;
//            auto c2 = ws.kci.compute_compression_for_node(m, 1);
//            if (std::isnan(c2) or c2<min_ci or c2>max_ci) continue;
//#pragma omp critical(inserting_hspnps)
//            {
//                //hl<<(n<m ? n:m)<<" "<<(n<m ? m:n)<<std::endl;
//                if (n < llabs(m)) hspnps.insert(std::make_pair(n, m));
//                else hspnps.insert(std::make_pair(llabs(m), (m>0 ? n:-n)));
//            }
//        }
//    }
//    return hspnps;
//}

//void LinkageUntangler::select_nodes_by_HSPNPs(uint64_t min_size, float min_ci, float max_ci) {
//
//    auto hspnps=get_HSPNPs(min_size,min_ci,max_ci);
//    sdglib::OutputLog() << "Selecting HSPNPs: " << hspnps.size() << " passed topology and CI" << std::endl;
//    for (auto p:hspnps) {
//        selected_nodes[llabs(p.first)] = true;
//        selected_nodes[llabs(p.second)] = true;
//    }
//
//}

void LinkageUntangler::select_linear_anchors(int min_links, int min_transitive_links) {
    for (auto n=1;n<dg.sdg.nodes.size();++n){
        if (dg.sdg.nodes[n].status == NodeStatus::Deleted) continue;
        bool anchor=true;
        auto fw_sorted=dg.fw_neighbours_by_distance(n,min_links);
        if (fw_sorted.size()>1) {
            for (auto i = 0; i < fw_sorted.size() - 1; ++i) {
                if (dg.link_count(-fw_sorted[i].second, fw_sorted[i + 1].second) < min_transitive_links) {
                    anchor = false;
                    break;
                }
            }
        }
        if (anchor) {
            auto bw_sorted = dg.fw_neighbours_by_distance(-n, min_links);
            if (bw_sorted.size()>1) {
                for (auto i = 0; i < bw_sorted.size() - 1; ++i) {
                    if (dg.link_count(-bw_sorted[i].second, bw_sorted[i + 1].second) < min_transitive_links) {
                        anchor = false;
                        break;
                    }
                }
            }
            if (bw_sorted.empty() and fw_sorted.empty()) anchor=false;
        }
        if (anchor) {
                selected_nodes[n] = true;
        }
    }
}

DistanceGraph LinkageUntangler::make_nextselected_linkage(int min_links) {
    DistanceGraph ldg(dg.sdg);
    for (auto n=1;n<dg.sdg.nodes.size();++n) {
        if (dg.sdg.nodes[n].status == NodeStatus::Deleted or !selected_nodes[n]) continue;
        auto fw_sorted = dg.fw_neighbours_by_distance(n, min_links);
        for (auto &fwl:fw_sorted) {
            if (llabs(fwl.second)!=n and selected_nodes[llabs(fwl.second)]) {
                if (!ldg.are_connected(-n, fwl.second)) {
                    ldg.add_link(-n, fwl.second, fwl.first);
                }
                break;
            }
        }
        auto bw_sorted = dg.fw_neighbours_by_distance(-n, min_links);
        for (auto &bwl:bw_sorted) {
            if (llabs(bwl.second)!=n and selected_nodes[llabs(bwl.second)]) {
                if (!ldg.are_connected(n, bwl.second)) {
                    ldg.add_link(n, bwl.second, bwl.first);
                }
                break;
            }
        }
    }
    return ldg;
}

//DistanceGraph LinkageUntangler::filter_linkage_to_hspnp_duos(uint64_t min_size, float min_ci, float max_ci,
//                                                              const DistanceGraph &ldg_old) {
//    std::unordered_map<sgNodeID_t,sgNodeID_t> node_to_parallel;
//    //1- get all hspnps -> create a map of parallels
//    DistanceGraph ldg_new(ws.sdg);
//    auto hspnps=get_HSPNPs(min_size,min_ci,max_ci);
//    for (auto h:hspnps) {
//        node_to_parallel[h.first]=h.second;
//        node_to_parallel[-h.first]=-h.second;
//        node_to_parallel[h.second]=h.first;
//        node_to_parallel[-h.second]=-h.first;
//    }
//    //2- hspnp -> look for links in one direction from one of the nodes, and same direction for the other
//    for (auto h:hspnps){
//        auto hr=h;
//        hr.first=-hr.first;
//        hr.second=-hr.second;
//        for (auto hspnp:{h,hr}) {
//            auto n1fs = ldg_old.get_fw_links(hspnp.first);
//            auto n2fs = ldg_old.get_fw_links(hspnp.second);
//            for (auto n1f:n1fs) {
//                for (auto n2f:n2fs) {
//                    if (node_to_parallel.count(n1f.dest) and node_to_parallel[n1f.dest] == n2f.dest) {
//                        // if links are to parts of the same node -> introduce linkage on newldg.
//                        ldg_new.add_link(-hspnp.first, n1f.dest, 0);
//                        ldg_new.add_link(-hspnp.second, n2f.dest, 0);
//                    }
//                }
//            }
//        }
//    }
//    return ldg_new;
//
//}

//void LinkageUntangler::expand_trivial_repeats(const DistanceGraph & ldg) {
//    uint64_t aa=0,ab=0,unsolved=0;
//    for (auto n=1;n<ws.sdg.nodes.size();++n) {
//        if (ws.sdg.nodes[n].status == NodeStatus::Deleted) continue;
//        //check node is 2-2
//        auto bwl=ws.sdg.get_bw_links(n);
//        if (bwl.size()!=2) continue;
//        auto fwl=ws.sdg.get_fw_links(n);
//        if (fwl.size()!=2) continue;
//        auto p1=-bwl[0].dest;
//        auto p2=-bwl[1].dest;
//        auto n1=fwl[0].dest;
//        auto n2=fwl[1].dest;
//        //check bw nodes have only one fw, is one of the fws and not the same
//        auto p1ll=ldg.get_fw_links(p1);
//        if (p1ll.size()!=1) continue;
//        auto p2ll=ldg.get_fw_links(p2);
//        if (p2ll.size()!=1) continue;
//        if (p1ll[0].dest==n1 and p2ll[0].dest==n2){
//            ws.sdg.expand_node(n,{{-p1},{-p2}},{{n1},{n2}});
//            ++aa;
//        }
//        else if (p2ll[0].dest==n1 and p1ll[0].dest==n2) {
//            ws.sdg.expand_node(n,{{-p1},{-p2}},{{n2},{n1}});
//            ++ab;
//        }
//        else ++unsolved;
//    }
//    sdglib::OutputLog()<<"Repeat expansion: AA:"<<aa<<"  AB:"<<ab<<"  Unsolved:"<<ab<<std::endl;
//}

//void LinkageUntangler::expand_linear_regions(const DistanceGraph & ldg) {
//    sdglib::OutputLog()<<"Starting linear region expansion..."<<std::endl;
//    //sdglib::OutputLog()<<"Looking for \"lines\"..."<<std::endl;
//    auto lines=ldg.get_all_lines(2);
//    sdglib::OutputLog()<<"Creating tag sets for "<<lines.size()<<" linear regions"<<std::endl;
//    //sdglib::OutputLog()<<"TODO: now use tags and LMPs to find paths between elements in the line"<<std::endl;
//    //sdglib::OutputLog()<<"USING ONLY 10 lines as a test"<<std::endl;
//    //lines.resize(10);
//    //---------------------------------Step 1: get tagsets for lines.
//    std::vector<std::set<bsg10xTag>> linetagsets;
//    linetagsets.reserve(lines.size());
//    BufferedTagKmerizer btk(ws.linked_reads_datastores[0],31,100000,1000);
//    for (auto l:lines){
//        //sdglib::OutputLog()<<"Analising line: ";
//        //for (auto &ln:l) std::cout<<"seq"<<llabs(ln)<<", ";
//        //for (auto &ln:l) std::cout<<ln<<" ";
//        //std::cout<<std::endl;
//        std::map<bsg10xTag ,std::pair<uint32_t , uint32_t >> tagcounts; //tag -> nodes, reads
//        for (auto &ln:l) {
//            std::map<bsg10xTag ,uint32_t> ntagcounts;
//            for (auto rm:ws.linked_reads_datastores[0].mapper.reads_in_node[llabs(ln)]){
//                auto tag=ws.linked_reads_datastores[0].get_read_tag(rm.read_id);
//                ++ntagcounts[tag];
//            }
//            for (auto ntc:ntagcounts) {
//                ++tagcounts[ntc.first].first;
//                tagcounts[ntc.first].second+=ntc.second;
//            }
//        }
//        std::map<bsg10xTag ,std::pair<uint32_t , uint32_t >> tagtotals;
//        std::set<bsg10xTag> lineTagSet;
//        for (auto tc:tagcounts) {
//            auto tag=tc.first;
//            auto reads=ws.linked_reads_datastores[0].get_tag_reads(tc.first);
//            std::set<sgNodeID_t> nodes;
//            for (auto r:reads) nodes.insert(ws.linked_reads_datastores[0].mapper.read_to_node[r]);
//            tagtotals[tag].first=nodes.size()-nodes.count(0);
//            tagtotals[tag].second=reads.size();
//            if (tc.second.first>1 and reads.size()<3000) lineTagSet.insert(tc.first);
//        }
//        linetagsets.push_back(lineTagSet);
//        if (linetagsets.size()%100==0) std::cout<<"."<<std::flush;
//    }
//    std::cout<<std::endl;
//    sdglib::OutputLog()<<"Creating path collections to be evaluated for "<<lines.size()<<" linear regions"<<std::endl;
//    std::vector<std::vector<std::vector<SequenceGraphPath>>> alternatives;
//    uint64_t total_paths=0,found=0,evaluated=0;
//    alternatives.reserve(lines.size());
//    for (auto l:lines) {
//        alternatives.emplace_back();
//        for (auto i = 0; i < l.size() - 1; ++i) {
//            evaluated++;
//            auto from = l[i];
//            auto to = l[i + 1];
//            auto paths = ws.sdg.find_all_paths_between(from, to, 400000, 20);
//            if (paths.size()>0) found++;
//            alternatives.back().emplace_back(paths);
//            total_paths+=paths.size();
//            //sdglib::OutputLog() << paths.size() << " paths to go from " << from << " to " << to << std::endl;
//        }
//        if (alternatives.size()%100==0) std::cout<<"."<<std::flush;
//    }
//    std::cout<<std::endl;
//    sdglib::OutputLog()<<"Junctions with possible paths: "<<found<<" / "<<evaluated<<std::endl;
//    sdglib::OutputLog()<<"Total paths to evaluate: "<<total_paths<<std::endl;
//    //Now use a function that only counts coverage on the set of kmers from all paths collections for each line
//    //kmer_coverage_in_tagreads(&std::map<kmer, coverage> (init at 0), std::set<tag> linetagset);
//
//    std::cout << "creating and populating the maps as of now" << std::endl;
//    std::vector<std::unordered_map<uint64_t, uint32_t>> linekmercoverages;
//    linekmercoverages.resize(lines.size());
//#pragma omp parallel
//    {
//        KmerMapCounter km_count(31);
//        KmerMapCreator km_create(31);
//        ReadSequenceBuffer blrsg(ws.linked_reads_datastores[0], 200000, 1000);
//        std::unordered_map<uint64_t, uint32_t> kmercoverages;
//        uint64_t done=0;
//#pragma omp for schedule(static, 100)
//        for (auto i = 0; i < lines.size(); ++i) {
//            //map with all kmers of paths to be evaluated
//            kmercoverages.clear();
////            size_t t=0;
////            for (auto &alts:alternatives[i]) {
////                for (auto &a:alts) {
////                    for (auto n:a.nodes) t += ws.sdg.nodes[llabs(n)].sequence.size();
////                }
////            }
////            kmercoverages.reserve(t);
//            for (auto &alts:alternatives[i]) {
//                for (auto &a:alts) {
//                    for (auto n:a.nodes) {
//                        km_create.create_all_kmers(ws.sdg.nodes[llabs(n)].sequence.c_str(), kmercoverages);
//                    }
//                }
//            }
//            for (auto &t:linetagsets[i]) {
//                for (auto rid:ws.linked_reads_datastores[0].get_tag_reads(t)) {
//                    km_count.count_all_kmers(blrsg.get_read_sequence(rid), kmercoverages);
//                }
//            }
//#pragma omp critical
//            linekmercoverages[i]=kmercoverages;
//            ++done;
//            if (done % 100 == 0) std::cout << "." << std::flush;
//            //count from the tag's reads
//            //btk.get_tag_kmers()
//        }
//    }
//    std::cout<<"DONE"<<std::endl;
//    sdglib::OutputLog()<<"evaluating alternative paths between each pair of adjacent nodes"<<std::endl;
//    KmerVectorCreator kvc(31);
//    uint64_t solved=0,none_covered=0,too_many_covered=0,no_paths=0;
//    GraphEditor ged(ws);
//    std::vector<SequenceGraphPath> sols;
//    for (auto i=0;i<lines.size();++i){
//        for (auto ia=0;ia<alternatives[i].size();++ia){
//            int best=-1;
//            bool too_many=false;
//            for (auto j=0;j<alternatives[i][ia].size();++j){
//                uint64_t missed=0;
//                for (auto n:alternatives[i][ia][j].nodes) {
//                    for (auto x:kvc.count_all_kmers(ws.sdg.nodes[llabs(n)].sequence.c_str())) {
//                        if (linekmercoverages[i][x] < 8) ++missed;//TODO: maybe ask for more than 1 read coverage?
//                    }
//                }
//                if (missed==0){
//                    if (best==-1) {
//                        best = j;
//                    }
//                    else {
//                        too_many=true;
//                        best = -1;
//                        break;
//                    }
//                }
//            }
//            //std::cout<<"Solution for line "<<i<<" jump #"<<ia<<": "<<best<<" / "<<alternatives[i][ia].size() <<std::endl;
//            if (best!=-1){
//                for (auto n:alternatives[i][ia][best].nodes) {
//                    if (selected_nodes[llabs(n)]){
//                        best=-1;
//                        break;
//                    }
//                }
//            }
//            if (best==-1) {
//                if (alternatives[i][ia].empty()) ++no_paths;
//                else if (too_many) ++too_many_covered;
//                else ++none_covered;
//            }
//            else {
//                ++solved;
//                sols.emplace_back(ws.sdg);
//                sols.back().nodes.emplace_back(lines[i][ia]);
//                for (auto n:alternatives[i][ia][best].nodes) sols.back().nodes.emplace_back(n);
//                sols.back().nodes.emplace_back(lines[i][ia+1]);
//            }
//        }
//    }
//    std::cout<<"Solved: "<<solved<<"  Too many covered paths: "<<too_many_covered<<"  No covered paths: "<<none_covered<<"  No paths found: "<<no_paths<<std::endl;
//    sdglib::OutputLog()<<"Applying solutions in the graph"<<std::endl;
//    uint64_t applied=0;
//    for (auto s:sols) {
//        if (ged.detach_path(s)) ++applied;
//    }
//    sdglib::OutputLog()<<applied<<" solutions applied"<<std::endl;
//}

//void LinkageUntangler::linear_regions_tag_local_assembly(const DistanceGraph & ldg, uint8_t k, int min_cvg, int max_lines, uint64_t min_nodes, uint64_t min_total_size, bool count_tag_cvg){
//    sdglib::OutputLog()<<"Starting linear region tag local assemblies..."<<std::endl;
//    auto lines=ldg.get_all_lines(min_nodes, min_total_size);
//    if (max_lines>0) {
//        sdglib::OutputLog()<<"USING ONLY "<<max_lines<< " lines out of "<<lines.size()<<" as a test"<<std::endl;
//        lines.resize(max_lines);
//        //to filter to specific node, keep only the line that has that node:
//        //auto old_lines=lines;
//        //lines.clear();
//        //for (auto l:lines) for (auto n:l) if (llabs(n)==TARGET_NODE_ID) lines.push_back(l);
//    }
//
//    sdglib::OutputLog()<<"Creating tag sets for "<<lines.size()<<" linear regions"<<std::endl;
//    //---------------------------------Step 1: get tagsets for lines.
//    std::vector<std::set<bsg10xTag>> linetagsets;
//    linetagsets.reserve(lines.size());
//    for (auto l:lines){
//        std::map<bsg10xTag ,std::pair<uint32_t , uint32_t >> tagcounts; //tag -> nodes, reads
//        for (auto &ln:l) {
//            std::map<bsg10xTag ,uint32_t> ntagcounts;
//            for (auto rm:ws.linked_reads_datastores[0].mapper.reads_in_node[llabs(ln)]){
//                auto tag=ws.linked_reads_datastores[0].get_read_tag(rm.read_id);
//                ++ntagcounts[tag];
//            }
//            for (auto ntc:ntagcounts) {
//                ++tagcounts[ntc.first].first;
//                tagcounts[ntc.first].second+=ntc.second;
//            }
//        }
//        std::map<bsg10xTag ,std::pair<uint32_t , uint32_t >> tagtotals;
//        std::set<bsg10xTag> lineTagSet;
//        for (auto tc:tagcounts) {
//            auto tag=tc.first;
//            auto reads=ws.linked_reads_datastores[0].get_tag_reads(tc.first);
//            std::set<sgNodeID_t> nodes;
//            for (auto r:reads) nodes.insert(ws.linked_reads_datastores[0].mapper.read_to_node[r]);
//            tagtotals[tag].first=nodes.size()-nodes.count(0);
//            tagtotals[tag].second=reads.size();
//            if (tc.second.first>1 and reads.size()<3000) lineTagSet.insert(tc.first);
//        }
//        linetagsets.push_back(lineTagSet);
//        if (linetagsets.size()%100==0) std::cout<<"."<<std::flush;
//    }
//    std::cout<<std::endl;
//
//    //-----Step 2: local assemblies
//    std::vector<SequenceGraphPath> sols;
//    std::atomic<uint64_t> found_transitions(0),not_found_transitions(0);
//    sdglib::OutputLog()<<"Performing local assembly for "<<lines.size()<<" linear regions"<<std::endl;
//    std::vector<std::vector<std::string>> local_unitigs;
//    local_unitigs.resize(lines.size());
//#pragma omp parallel
//    {
//        ReadSequenceBuffer blrsg(ws.linked_reads_datastores[0], 200000, 1000);
//        std::vector<SequenceGraphPath> tsols;
//        uint64_t donelines = 0;
//#pragma omp for schedule(dynamic, 1)
//        for (auto i = 0; i < lines.size(); ++i) {
//            auto ltkmers128 = ws.linked_reads_datastores[0].get_tags_kmers128(k, min_cvg, linetagsets[i], blrsg,
//                                                                             count_tag_cvg);
//            //std::cout << "creating DBG for line #" << i << std::endl;
//            WorkSpace pws;
//            SequenceDistanceGraph dbg(pws);
//            GraphMaker gm(dbg);
//            gm.new_graph_from_kmerset_trivial128(ltkmers128, k);
//            //dbg.write_to_gfa("local_dbg_" + std::to_string(i) + ".gfa");
//            //gruesome tip clipping:
//            //std::cout << "Starting gruesome tip clipping" << std::endl;
//            std::set<sgNodeID_t> to_delete;
//            for (sgNodeID_t n = 1; n < dbg.nodes.size(); ++n) {
//                if (dbg.nodes[n].status == NodeStatus::Deleted) continue;
//                if (dbg.nodes[n].sequence.size() > 200) continue;
//                //std::cout<<"Evaluating seq"<<n<<": ";
//                auto fwl = dbg.get_fw_links(n);
//                auto bwl = dbg.get_bw_links(n);
//                //std::cout<<" fwl: "<<fwl.size()<<"  bwl: "<<bwl.size();
//                if (fwl.size() == 1 and bwl.size() == 0) {
//                    //std::cout<<"  bwl for "<<fwl[0].dest<<": "<<dbg.get_bw_links(fwl[0].dest).size();
//                    if (dbg.get_bw_links(fwl[0].dest).size() == 2) {
//                        to_delete.insert(n);
//                        //std::cout<<" D"<<std::endl;
//                    }
//                }
//                if (fwl.size() == 0 and bwl.size() == 1) {
//                    //std::cout<<"  fwl for "<<-bwl[0].dest<<": "<<dbg.get_fw_links(-bwl[0].dest).size();
//                    if (dbg.get_fw_links(-bwl[0].dest).size() == 2) {
//                        to_delete.insert(n);
//                        //std::cout<<" D"<<std::endl;
//                    }
//                }
//                if (fwl.size() == 0 and bwl.size() == 0) to_delete.insert(n);
//                //std::cout<<std::endl;
//            }
//            //std::cout << "Nodes to delete: " << to_delete.size() << std::endl;
//            for (auto n:to_delete) dbg.remove_node(n);
//            auto utc = dbg.join_all_unitigs();
//#pragma omp critical
//            {
//
//
//                for (auto n = 0; n < dbg.nodes.size(); ++n) {
//                    if (dbg.nodes[n].sequence.size() > 2000) {
//                        //patch_unitigs << ">local_dbg_" << i << "_node_" << n << std::endl << dbg.nodes[n].sequence
//                        //              << std::endl;
//                        local_unitigs[i].emplace_back(dbg.nodes[n].sequence);
//                    }
//                }
//            }
//        }
//    }
//    sdglib::OutputLog()<<"Local assemblies done!"<<std::endl;
//    //---------- Step 3: patching graph
//    sdglib::OutputLog()<<"Expanding local unitigs to their RC"<<std::endl;
//    //first expand local unitigs into direct and rc
//
//    for (auto &lu:local_unitigs) {
//        auto luc=lu.size();
//        lu.reserve(2*luc);
//        for (auto i = 0; i < luc; ++i) {
//            Node n(lu[i]);
//            n.make_rc();
//            lu.emplace_back(n.sequence);
//        }
//    }
//    sdglib::OutputLog()<<"Looking for transtitions to patch"<<std::endl;
//    std::vector<std::pair<std::pair<sgNodeID_t , sgNodeID_t >, std::string>> patches;
//    GraphEditor ge(ws);
//#pragma omp parallel for
//        for (auto i = 0; i < lines.size(); ++i) {
//            //sdglib::OutputLog()<<"Line #"<<i<<std::endl;
//            for (auto li = 0; li < lines[i].size() - 1; ++li) {
//                auto n1 = ws.sdg.nodes[llabs(lines[i][li])];
//                auto n2 = ws.sdg.nodes[llabs(lines[i][li + 1])];
//                const size_t ENDS_SIZE = 1000;
//                if (n1.sequence.size() > ENDS_SIZE)
//                    n1.sequence = n1.sequence.substr(n1.sequence.size() - ENDS_SIZE - 1, ENDS_SIZE);
//                if (n2.sequence.size() > ENDS_SIZE) n2.sequence.resize(ENDS_SIZE);
//                if (lines[i][li] < 0) n1.make_rc();
//                if (lines[i][li + 1] < 0) n2.make_rc();
//                std::vector<std::string> matches;
//                for (auto &unitig:local_unitigs[i]) {
//                    auto n1pos = unitig.find(n1.sequence);
//                    auto n2pos = unitig.find(n2.sequence);
//                    if (n1pos < unitig.size() and n2pos < unitig.size()) {
//                        //std::cout << lines[i][li] << " and " << lines[i][li + 1] << " found on unitig " << n
//                        //          << std::endl;
//                        matches.emplace_back(unitig.substr(n1pos, n2pos + 2 * ENDS_SIZE - n1pos));
//                    }
//                }
//                //TODO: collapse unitigs that are equivalent.
//                if (matches.size() == 1) {
//#pragma omp critical
//                    //patches.emplace_back(std::make_pair(lines[i][li],lines[i][li+1]),matches[0]);
////                    //std::cout<<" Patching between "<<lines[i][li]<<" and "<<lines[i][li+1]<<std::endl;
////                    auto prevn = ws.sdg.nodes.size();
//                    auto patch_code = ge.patch_between(lines[i][li], lines[i][li + 1], matches[0]);
//////                    patchfile << ">patch_" << lines[i][li] << "_" << lines[i][li + 1] << "_"
//////                              << (ws.sg.nodes.size() > prevn ? "APPLIED_" : "FAILED_")
//////                              << patch_code << std::endl << matches[0] << std::endl;
////                    //std::cout<<" Patched!!!"<<std::endl;
////                    if (ws.sdg.nodes.size() > prevn) ++patched;
////                    else ++not_patched;
//                }
//
//            }
//        }
//    sdglib::OutputLog()<<"Joining unitigs"<<std::endl;
//    auto ujc=ws.sdg.join_all_unitigs();
//    sdglib::OutputLog()<<"Unitigs joined after patching: "<<ujc<<std::endl;
//
//}

//void LinkageUntangler::expand_linear_regions_skating(const DistanceGraph & ldg, int max_lines) {
//    sdglib::OutputLog()<<"Starting linear region consolidation via skating with line tag collection..."<<std::endl;
//    auto lines=ldg.get_all_lines(2);
//    if (max_lines>0) {
//        sdglib::OutputLog()<<"USING ONLY "<<max_lines<< " lines as a test"<<std::endl;
//        lines.resize(max_lines);
//    }
//
//    sdglib::OutputLog()<<"Creating tag sets for "<<lines.size()<<" linear regions"<<std::endl;
//    //---------------------------------Step 1: get tagsets for lines.
//    std::vector<std::set<bsg10xTag>> linetagsets;
//    linetagsets.reserve(lines.size());
//    BufferedTagKmerizer btk(ws.linked_reads_datastores[0],31,100000,1000);
//    for (auto l:lines){
//        //sdglib::OutputLog()<<"Analising line: ";
//        //for (auto &ln:l) std::cout<<"seq"<<llabs(ln)<<", ";
//        //for (auto &ln:l) std::cout<<ln<<" ";
//        //std::cout<<std::endl;
//        std::map<bsg10xTag ,std::pair<uint32_t , uint32_t >> tagcounts; //tag -> nodes, reads
//        for (auto &ln:l) {
//            std::map<bsg10xTag ,uint32_t> ntagcounts;
//            for (auto rm:ws.linked_reads_datastores[0].mapper.reads_in_node[llabs(ln)]){
//                auto tag=ws.linked_reads_datastores[0].get_read_tag(rm.read_id);
//                ++ntagcounts[tag];
//            }
//            for (auto ntc:ntagcounts) {
//                ++tagcounts[ntc.first].first;
//                tagcounts[ntc.first].second+=ntc.second;
//            }
//        }
//        std::map<bsg10xTag ,std::pair<uint32_t , uint32_t >> tagtotals;
//        std::set<bsg10xTag> lineTagSet;
//        for (auto tc:tagcounts) {
//            auto tag=tc.first;
//            auto reads=ws.linked_reads_datastores[0].get_tag_reads(tc.first);
//            std::set<sgNodeID_t> nodes;
//            for (auto r:reads) nodes.insert(ws.linked_reads_datastores[0].mapper.read_to_node[r]);
//            tagtotals[tag].first=nodes.size()-nodes.count(0);
//            tagtotals[tag].second=reads.size();
//            if (tc.second.first>1 and reads.size()<3000) lineTagSet.insert(tc.first);
//        }
//        linetagsets.push_back(lineTagSet);
//        if (linetagsets.size()%100==0) std::cout<<"."<<std::flush;
//    }
//    std::cout<<std::endl;
//    uint64_t jc=0;
//    for (auto &l:lines) jc+=l.size()-1;
//    sdglib::OutputLog()<<"Skating across "<<jc<<" junctions in "<<lines.size()<<" linear regions"<<std::endl;
//
//    std::vector<SequenceGraphPath> sols;
//#pragma omp parallel
//    {
//        ReadSequenceBuffer blrsg(ws.linked_reads_datastores[0], 200000, 1000);
//        std::vector<SequenceGraphPath> tsols;
//        uint64_t donelines=0;
//#pragma omp for schedule(dynamic,1)
//        for (auto i=0; i<lines.size(); ++i){
//            //std::cout<<"Creating kmer set for line"<<i<<" from tags"<<std::endl;
//            auto ltkmers=ws.linked_reads_datastores[0].get_tags_kmers(31,3,linetagsets[i],blrsg);
//            UncoveredKmerCounter ukc(31,ltkmers);
//            //std::cout<<"Evaluating paths for "<<lines[i].size()-1<<" junctions"<<std::endl;
//            for (auto j=0;j<lines[i].size()-1;++j){
//                auto from=lines[i][j];
//                auto to=lines[i][j+1];
//                //std::cout<<std::endl<<std::endl<<"Junction #"<<j+1<<" from "<<from<<" to "<<to<<std::endl;
//                std::vector<std::vector<sgNodeID_t>> skated_paths;
//                skated_paths.push_back({from});
//                int max_nodes=50;
//                while (--max_nodes and not skated_paths.empty()){
//                    //std::cout<<std::endl<<"expansion round starting with "<<skated_paths.size()<<" paths "<<std::endl;
//                    auto old_skated=skated_paths;
//                    skated_paths.clear();
//                    bool loop=false,crosstalk=false;
//                    for (auto p:old_skated) {
//                        if (p.back()==to) {
//                            skated_paths.push_back(p);
//                            continue;
//                        }
//                        //std::cout<<" expanding fw from node "<<p.back()<<std::endl;
//                        for (auto fwl:ws.sdg.get_fw_links(p.back())) {
//                            //std::cout<<"  considering fwl to "<<fwl.dest<<std::endl;
//                            if (std::count(p.begin(),p.end(),fwl.dest)>0 or std::count(p.begin(),p.end(),-fwl.dest)>0){
//                                loop=true;
//                                //std::cout<<"loop detected, aborting junction analysis"<<std::endl;
//                                break;
//                            }
//
//                            auto u=ukc.count_uncovered(ws.sdg.nodes[llabs(fwl.dest)].sequence.c_str());
//                            //std::cout<<"  Uncovered kmers in "<<fwl.dest<<" ("<<ws.sdg.nodes[llabs(fwl.dest)].sequence.size()<<" bp): "
//                            //                                                                                                <<u<<std::endl;
//                            if ( u == 0) {
//                                //check for a path that reaches a selected node that is not connected here
//                                if (selected_nodes[llabs(fwl.dest)] and fwl.dest!=to) {
//                                    crosstalk=true;
//                                    break;
//                                }
//                                //std::cout<<"  path can continue in node"<<fwl.dest<<std::endl;
//                                skated_paths.push_back(p);
//                                skated_paths.back().push_back(fwl.dest);
//                            }
//                        }
//                    }
//                    if (loop or crosstalk) {
//                        skated_paths.clear();
//                        break;
//                    }
//                }
//                uint64_t complete=0,incomplete=0;
//                for (auto p:skated_paths) {
//                    if (p.back()==to) ++complete;
//                    else ++incomplete;
//                }
//                if (complete==1 and incomplete==0) tsols.emplace_back(SequenceGraphPath(ws.sdg,skated_paths[0]));
//                //std::cout<<"Skating line #"<<i+1<<" junction #"<<j+1<<" produced "<<complete<<" complete paths and "<<incomplete<<" possibly incomplete paths"<<std::endl;
//            }
//            if (++donelines%100==0) std::cout<<"."<<std::flush;
//        }
//#pragma omp critical
//        sols.insert(sols.end(),tsols.begin(),tsols.end());
//    }
//    sdglib::OutputLog()<<"Applying "<<sols.size()<<" solutions in the graph"<<std::endl;
//    GraphEditor ged(ws);
//    uint64_t applied=0;
//    for (auto s:sols) {
//        if (ged.detach_path(s)) ++applied;
//    }
//    sdglib::OutputLog()<<applied<<" solutions applied"<<std::endl;
//}

//void LinkageUntangler::fill_linkage_line(std::vector<sgNodeID_t> nodes) {
//    std::cout<<"Filling linkage for line:";
//    for (auto n:nodes) std::cout<<" "<<n;
//    std::cout<<std::endl;
//    std::cout<<"Creating a set with all the possibly local reads"<<std::endl;
//    //TODO: create the set, with both 10x and LMP/PE reads
//    std::map<bsg10xTag ,std::pair<uint32_t , uint32_t >> tagcounts; //tag -> nodes, reads
//    for (auto &ln:nodes) {
//        std::map<bsg10xTag ,uint32_t> ntagcounts;
//        for (auto rm:ws.linked_reads_datastores[0].mapper.reads_in_node[llabs(ln)]){
//            auto tag=ws.linked_reads_datastores[0].get_read_tag(rm.read_id);
//            if (tag==0) continue;
//            ++ntagcounts[tag];
//        }
//        for (auto ntc:ntagcounts) {
//            ++tagcounts[ntc.first].first;
//            tagcounts[ntc.first].second+=ntc.second;
//        }
//    }
//    std::map<bsg10xTag ,std::pair<uint32_t , uint32_t >> tagtotals;
//    std::set<bsg10xTag> lineTagSet;
//    uint64_t total_reads=0;
//    for (auto tc:tagcounts) {
//        auto tag=tc.first;
//        auto reads=ws.linked_reads_datastores[0].get_tag_reads(tc.first);
//        total_reads+=reads.size();
//        std::set<sgNodeID_t> nodes;
//        for (auto r:reads) nodes.insert(ws.linked_reads_datastores[0].mapper.read_to_node[r]);
//        tagtotals[tag].first=nodes.size()-nodes.count(0);
//        tagtotals[tag].second=reads.size();
//        if (tc.second.first>1 and reads.size()<3000) lineTagSet.insert(tc.first);
//    }
//    std::cout<<"Local tag reads: "<<total_reads<<std::endl;
//    std::cout<<"Creating an uncleaned DBG"<<std::endl;
//    ReadSequenceBuffer blrsg(ws.linked_reads_datastores[0], 200000, 1000);
//    auto ltkmers128 = ws.linked_reads_datastores[0].get_tags_kmers128(63, 3, lineTagSet, blrsg, true);
//    //std::cout << "creating DBG for line #" << i << std::endl;
//    WorkSpace pws;
//    SequenceDistanceGraph dbg(pws);
//    GraphMaker gm(dbg);
//    gm.new_graph_from_kmerset_trivial128(ltkmers128, 63);
//    std::ofstream anchf("local_dbg_" + std::to_string(nodes[0]) + "_anchors.fasta");
//    for (auto n:nodes){
//        anchf<<">seq"<<llabs(n)<<std::endl;
//        anchf<<ws.sdg.nodes[llabs(n)].sequence<<std::endl;
//
//    }
//    dbg.write_to_gfa1("local_dbg_"+std::to_string(nodes[0])+"_uncleaned.gfa");
//    gm.tip_clipping(200);
//    gm.remove_small_unconnected(500);
//    dbg.write_to_gfa1("local_dbg_"+std::to_string(nodes[0])+".gfa");
//
//    std::cout<<"Analising junctions, one by one"<<std::endl;
//    for (auto i=0;i<nodes.size()-1;++i){
//        std::cout<<"Tring to joing "<<nodes[i]<<" (-) -> (+) "<<nodes[i+1]<<std::endl;
//    }
//    /*
//
//    linetagsets.push_back(lineTagSet);
//    if (linetagsets.size()%100==0) std::cout<<"."<<std::flush;*/
//}
