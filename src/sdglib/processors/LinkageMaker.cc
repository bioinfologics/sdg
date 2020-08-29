//
// Created by Bernardo Clavijo (EI) on 2019-07-17.
//

#include <atomic>
#include "GraphMaker.hpp"
#include "LinkageUntangler.hpp"
#include "LinkageMaker.hpp"
#include <sdglib/mappers/LongReadsMapper.hpp>
#include <sdglib/mappers/LinkedReadsMapper.hpp>

LongReadHaplotypeMappingsFilter::LongReadHaplotypeMappingsFilter (const LongReadsMapper & _lorm, const LinkedReadsMapper & _lirm):
        lorm(_lorm),lirm(_lirm){
    lrbsgp=new ReadSequenceBuffer(lorm.datastore);
}

void LongReadHaplotypeMappingsFilter::set_read(uint64_t _read_id) {
    read_id=_read_id;
    read_seq=std::string(lrbsgp->get_read_sequence(_read_id));
    mappings = lorm.get_raw_mappings_from_read(read_id);
    //TODO: filter mappings?
    nodeset.clear();
    for (auto &m:mappings) {
        auto n=llabs(m.node);
        bool dup=false;
        for (auto &x:nodeset) if (x==n) { dup=true; break; }
        if (!dup) nodeset.emplace_back(n);
    }
    std::sort(nodeset.begin(),nodeset.end());
    haplotype_scores.clear();
    rkmers.clear();
}

void LongReadHaplotypeMappingsFilter::generate_haplotypes_from_linkedreads(float min_tn) {
    //2) create the neighbour sets, uniq to not evaluate same set many times over.
    haplotype_scores.clear();
    std::vector<std::vector<sgNodeID_t>> all_nsets;
    for (auto n:nodeset){
        all_nsets.emplace_back();
        if (n>=lirm.tag_neighbours.size() or lirm.tag_neighbours[n].empty()) continue;
        for (auto m:nodeset){
            for (auto &score:lirm.tag_neighbours[n]){
                if (score.node==m and score.score>min_tn) {
                    all_nsets.back().push_back(m);
                    break;
                }
            }
        }
    }
    std::sort(all_nsets.begin(),all_nsets.end());
    for (const auto &ns:all_nsets) {
        if (ns.empty()) continue;
        bool dup=false;
        for (auto &hs:haplotype_scores) if (ns==hs.haplotype_nodes) {dup=true;break;}
        if (dup) continue;
        haplotype_scores.emplace_back();
        haplotype_scores.back().haplotype_nodes=ns;
        haplotype_scores.back().score=0;
    }
}

void LongReadHaplotypeMappingsFilter::score_coverage(float weight) {
    for (auto &hs:haplotype_scores){
        coverage.clear();
        coverage.resize(read_seq.size());
        auto &nodeset=hs.haplotype_nodes;
        uint64_t total_map_bp=0;
        for (auto m:mappings) if (std::count(nodeset.begin(),nodeset.end(),llabs(m.node))) total_map_bp+=m.qEnd-m.qStart;
        //this can be done faster by saving all starts and ends and keeping a rolling value;
        for (auto m:mappings) {
            if (std::count(nodeset.begin(),nodeset.end(),llabs(m.node))){
                for (auto i=m.qStart;i<=m.qEnd and i<read_seq.size();++i) ++coverage[i];
            }
        }
        int64_t cov1=std::count(coverage.begin(),coverage.end(),1);
        int64_t cov0=std::count(coverage.begin(),coverage.end(),0);
        int64_t covm=coverage.size()-cov0-cov1;
        float score=(cov1-4*covm )*1.0/coverage.size();
        hs.score+=weight*score;
    }
}

void LongReadHaplotypeMappingsFilter::score_window_winners(float weight, int k, int win_size, int win_step) {
    //kmerise read
    StreamKmerFactory skf(k);
    skf.produce_all_kmers(read_seq.c_str(),rkmers);
    std::unordered_map<uint64_t,int32_t> rkmerpos;
    rkmerpos.reserve(rkmers.size());
    int32_t rkp=0;
    for (auto &rk:rkmers){
        auto rkhit=rkmerpos.find(rk);
        if (rkhit!=rkmerpos.end()) rkmerpos[rk]=-1;
        else rkmerpos[rk]=rkp;
        ++rkp;
    };
    if (rkmers.size()<win_size) return;
    if (kmer_hits_by_node.size()<nodeset.size()) kmer_hits_by_node.resize(nodeset.size());
    //first create a list of used kmers
    uint64_t ni=0;
    for (auto &n:nodeset){
        nkmers.clear();
        skf.produce_all_kmers(lorm.sg.nodes[n].sequence.c_str(),nkmers);
        //std::sort(nkmers.begin(),nkmers.end());
        kmer_hits_by_node[ni].clear();
        kmer_hits_by_node[ni].resize(rkmers.size());
        for (auto &nk:nkmers) {
            auto rkhit=rkmerpos.find(nk);
            if (rkhit!=rkmerpos.end() and rkhit->second!=-1) kmer_hits_by_node[ni][rkhit->second]=true;
        }
        /*for (auto i=0;i<rkmers.size();++i){
            auto nklb=std::lower_bound(nkmers.begin(),nkmers.end(),rkmers[i]);
            if (nklb!=nkmers.end() and *nklb==rkmers[i]) kmer_hits_by_node[ni][i]=1;
        }*/
        ++ni;
    }
    std::map<sgNodeID_t,int> node_wins;
    int total_windows=0;
    int total_wins=0;
    for (int64_t start=0;start<rkmers.size()-win_size;start+=win_step,++total_windows){
        int win_node_idx=-1;
        int win_score=.2*win_size;
        for (auto ni=0;ni<nodeset.size();++ni){
            auto nscore=std::count(kmer_hits_by_node[ni].begin()+start,kmer_hits_by_node[ni].begin()+start+win_size,1);
            if (nscore>win_score) {
                win_node_idx=ni;
                win_score=nscore;
            }
        }
        if (win_node_idx!=-1) {
            ++total_wins;
            ++node_wins[nodeset[win_node_idx]];
        }
    }

    for (auto &hs:haplotype_scores){
        float score=0;
        for (auto &n:hs.haplotype_nodes) score+=node_wins[n];
        score/=total_windows;
        hs.score+=weight*score;
    }
}

void LongReadHaplotypeMappingsFilter::rank(uint64_t read_id, float coverage_weight, float winners_weight, float min_tn) {
    set_read(read_id);
    generate_haplotypes_from_linkedreads(min_tn);
    score_window_winners(winners_weight);
    score_coverage(coverage_weight);
}

void LinkageMaker::select_all() {
    selected_nodes.clear();
    selected_nodes.resize(dg.sdg.nodes.size(),true);
}

void LinkageMaker::deselect_all() {
    selected_nodes.clear();
    selected_nodes.resize(dg.sdg.nodes.size(),false);
}

void LinkageMaker::report_selection() {
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
void LinkageMaker::select_by_size( uint64_t min_size, uint64_t max_size) {
    for (auto n=1;n<dg.sdg.nodes.size();++n) {
        if (dg.sdg.nodes[n].status==NodeStatus::Deleted) continue;
        if (dg.sdg.nodes[n].sequence.size() >= min_size and
            (max_size==0 or dg.sdg.nodes[n].sequence.size() <= max_size))
            selected_nodes[n]=true;
    }
}

DistanceGraph LinkageMaker::make_topology_linkage(int radius) {
    DistanceGraph ldg(dg.sdg);
    for (auto m=1; m < dg.sdg.nodes.size(); ++m) {
        if (!selected_nodes[m]) continue;
        for (auto n:{m,-m}) {
            std::set<sgNodeID_t> reached, last = {n};
            for (auto i = 0; i < radius; ++i) {
                std::set<sgNodeID_t> new_last;
                for (auto l:last) {
                    for (auto fwl:dg.sdg.get_fw_links(l)) {
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

DistanceGraph LinkageMaker::make_paired_linkage(int min_reads) {
    DistanceGraph ldg(dg.sdg);
    std::map<std::pair<sgNodeID_t, sgNodeID_t>, uint64_t> lv;
    sdglib::OutputLog()<<"collecting link votes across all paired libraries"<<std::endl;
    //use all libraries collect votes on each link
    auto rmi=0;
    for (auto &prds:dg.sdg.ws.paired_reads_datastores) {
        for (auto i = 1; i < prds.mapper.read_to_node.size(); i += 2) {
            sgNodeID_t n1 = prds.mapper.read_to_node[i];
            sgNodeID_t n2 = prds.mapper.read_to_node[i + 1];
            if (n1 == 0 or n2 == 0 or n1 == n2 or !selected_nodes[n1] or !selected_nodes[n2] ) continue;
            if (prds.mapper.read_direction_in_node[i]) n1=-n1;
            if (prds.mapper.read_direction_in_node[i+1]) n2=-n2;
            if (llabs(n1) > llabs(n2)) std::swap(n1,n2);
            ++lv[std::make_pair(n1, n2)];
        }
        ++rmi;
    }

    sdglib::OutputLog()<<"adding links with "<<min_reads<<" votes"<<std::endl;
    //std::vector<std::vector<std::pair<sgNodeID_t ,uint64_t>>> nodelinks(ws.sdg.nodes.size());
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

DistanceGraph LinkageMaker::make_paired_linkage_pe(int min_reads) {
    DistanceGraph ldg(dg.sdg);
    std::map<std::pair<sgNodeID_t, sgNodeID_t>, uint64_t> lv;
    sdglib::OutputLog()<<"collecting link votes across all paired libraries"<<std::endl;
    //use all libraries collect votes on each link
    auto &prds= dg.sdg.ws.paired_reads_datastores[0];
    for (auto i = 1; i < prds.mapper.read_to_node.size(); i += 2) {
        sgNodeID_t n1 = prds.mapper.read_to_node[i];
        sgNodeID_t n2 = prds.mapper.read_to_node[i + 1];
        if (n1 == 0 or n2 == 0 or n1 == n2 or !selected_nodes[n1] or !selected_nodes[n2] ) continue;
        if (not prds.mapper.read_direction_in_node[i]) n1=-n1;
        if (not prds.mapper.read_direction_in_node[i+1]) n2=-n2;
        if (llabs(n1) > llabs(n2)) std::swap(n1,n2);
        ++lv[std::make_pair(n1, n2)];
    }


    sdglib::OutputLog()<<"adding links with "<<min_reads<<" votes"<<std::endl;
    //std::vector<std::vector<std::pair<sgNodeID_t ,uint64_t>>> nodelinks(ws.sdg.nodes.size());
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
void add_readkmer_nodes(std::vector<sgNodeID_t> & kmernodes, std::vector<std::pair<uint64_t,bool>> & readkmers, const std::unordered_map<uint64_t, graphStrandPos> & index, bool rev){
    //TODO allow for a minimum of kmers to count the hit?
    if (not rev) {
        for (auto rki=readkmers.begin();rki<readkmers.end();++rki) {
            auto kp=index.find(rki->first);
            if (kp==index.end()) continue; //kmer not found
            auto node=(rki->second ? -kp->second.node:kp->second.node);
            if (kmernodes.empty() or kmernodes.back()!=node) kmernodes.emplace_back(node);
        }
    }
    else {
        for (auto rki=readkmers.rbegin();rki<readkmers.rend();++rki) {
            auto kp=index.find(rki->first);
            if (kp==index.end()) continue; //kmer not found
            auto node=(rki->second ? kp->second.node:-kp->second.node);
            if (kmernodes.empty() or kmernodes.back()!=node) kmernodes.emplace_back(node);
        }
    }

}

std::map<std::pair<sgNodeID_t, sgNodeID_t>, uint64_t> LinkageMaker::shared_read_paths(int min_shared, std::vector<size_t> libraries, bool r1rev, bool r2rev ) {
    std::map<std::pair<sgNodeID_t, sgNodeID_t>, uint64_t> shared_paths;
    std::cout<<"Creating paired linkage by kmer"<<std::endl;
    auto ukindex=UniqueKmerIndex(dg.sdg);
    std::cout<<"Looking up reads"<<std::endl;
    std::vector<std::pair<sgNodeID_t ,sgNodeID_t >> nodeproximity;
    //This actually works like a paired-read-to-path
    for (auto lib:libraries) {
#pragma omp parallel
        {
            CStringKMerFactory cskf(31);
            std::vector<std::pair<sgNodeID_t ,sgNodeID_t >> nodeproximity_thread;
            std::vector<std::pair<uint64_t,bool>> read1kmers,read2kmers;
            std::vector<sgNodeID_t> kmernodes;
            ReadSequenceBuffer bprsg(dg.sdg.ws.paired_reads_datastores[lib], 1000000, dg.sdg.ws.paired_reads_datastores[lib].readsize * 2 + 2);
#pragma omp for
            for (auto rid = 1; rid < dg.sdg.ws.paired_reads_datastores[lib].size(); rid += 2) {
                //std::cout<<"analising reads "<<rid<<" and "<<rid+1<<std::endl;

                read1kmers.clear();
                read2kmers.clear();
                kmernodes.clear();

                cskf.create_kmers_direction(read1kmers, bprsg.get_read_sequence(rid));
                cskf.create_kmers_direction(read2kmers, bprsg.get_read_sequence(rid + 1));
                //first put the kmers from read 1 in there;
                add_readkmer_nodes(kmernodes, read1kmers, ukindex.getMap(), r1rev);
                //for (auto kn:kmernodes) std::cout<<" "<<kn; std::cout<<std::endl;
                add_readkmer_nodes(kmernodes, read2kmers, ukindex.getMap(), r2rev);
                //for (auto kn:kmernodes) std::cout<<" "<<kn; std::cout<<std::endl;
                //TODO: A bit of a more clever connection thingy (i.e. off to errors and back?)
                for (auto i = 0; i < kmernodes.size(); ++i) {
                    auto &kn1 = kmernodes[i];
                    for (auto j = i + 1; j < kmernodes.size(); ++j) {
                        auto &kn2 = kmernodes[j];
                        if (llabs(kn1) < llabs(kn2)) {
                            nodeproximity_thread.emplace_back(-kn1, kn2);
                        } else {
                            nodeproximity_thread.emplace_back(kn2, -kn1);
                        }
                    }
                }

            }
#pragma omp critical(collect_nodeproximity)
            {
               nodeproximity.insert(nodeproximity.end(),nodeproximity_thread.begin(),nodeproximity_thread.end());
               nodeproximity_thread.clear();
            }
        }
    }
    std::cout<<"Collecting proximity totals"<<std::endl;
    std::ofstream ptotf("proximity_detail.csv");
    if (not nodeproximity.empty()) {
        std::sort(nodeproximity.begin(), nodeproximity.end());
        std::pair<sgNodeID_t, sgNodeID_t> curr_pair = {0, 0};
        size_t curr_first = 0;
        size_t i;
        for (i = 0; i < nodeproximity.size(); ++i) {
            if (nodeproximity[i] != curr_pair) {
                auto c = i - curr_first;
                if (c >= min_shared) {
                    shared_paths[curr_pair]=c;
                    ptotf << curr_pair.first << ", " << curr_pair.second << ", " << c << std::endl;
                }
                curr_first = i;
                curr_pair = nodeproximity[i];
            }
        }
        auto c = i - curr_first;
        if (c >= min_shared) {
            shared_paths[curr_pair]=c;
            ptotf << curr_pair.first << ", " << curr_pair.second << ", " << c << std::endl;
        }
    }
    std::cout<<"Done!"<<std::endl;
    //TODO: create linkage
    return shared_paths;
}

std::vector<Link> LinkageMaker::mappings_to_multilinkage(const std::vector<LongReadMapping> &lorm_mappings, uint32_t read_size, int32_t unmapped_end) {
    std::vector<Link> linkage;
    std::unordered_map<sgNodeID_t,uint64_t> total_bp;
    std::vector<LongReadMapping> mfilt_total,mmergedfilt;
    //first, copy mappings, removing nodes whose total mapping adds up to less than 70% of the node unless first or last
    for (auto &m:lorm_mappings) total_bp[m.node] +=m.nEnd-m.nStart+1;
    if (read_size==0) for (auto &m:lorm_mappings) if (m.qEnd>read_size) read_size=m.qEnd;
    for (int i=0; i<lorm_mappings.size();++i){
        auto &m=lorm_mappings[i];
        if (not selected_nodes[llabs(m.node)]) continue;
        auto ns= dg.sdg.nodes[llabs(m.node)].sequence.size();
        if ( (i==0 and m.qStart<unmapped_end) or (i==lorm_mappings.size()-1 and m.qEnd+unmapped_end>=read_size) or total_bp[m.node]>=.7*ns) {
            //If a node has mode than one consecutive mapping, merge them.
            if (mfilt_total.size()>0 and mfilt_total.back().node==m.node and mfilt_total.back().nStart<m.nStart and mfilt_total.back().nEnd<m.nEnd and mfilt_total.back().nEnd<m.nStart+500) {
                mfilt_total.back().nEnd = m.nEnd;
                mfilt_total.back().qEnd = m.qEnd;
            }
            else mfilt_total.push_back(m);
        }
    }
    //Now remove all mappings that do not cover 80% of the node
    for (int i=0; i<mfilt_total.size();++i){
        auto &m=mfilt_total[i];
        auto ns= dg.sdg.nodes[llabs(m.node)].sequence.size();
        if ( (i==0 and m.nEnd>.9*ns and m.qStart<unmapped_end) or (i==mfilt_total.size()-1 and m.nStart<.1*ns and m.qEnd+unmapped_end>read_size) or m.nEnd-m.nStart+1>=.8*ns) {
            mmergedfilt.push_back(m);
        }
    }
    //now compute starts and ends for 0% and 100%
    std::vector<std::pair<sgNodeID_t, std::pair<int32_t, int32_t>>> node_ends;
    for (auto &m:mmergedfilt) {
        node_ends.emplace_back(m.node, std::make_pair(m.qStart-m.nStart, m.qEnd +
                                                                         dg.sdg.nodes[llabs(m.node)].sequence.size() - m.nEnd));
    }
    //for every nodeA:
    for (int nA=0;nA+1<node_ends.size();++nA) {
        //for every other nodeB fw:
        for (int nB=nA+1;nB<node_ends.size();++nB) {
            //link from -nodeA to +nodeB with dist start[nodeB]-end[nodeA]-1
            linkage.emplace_back(-node_ends[nA].first,node_ends[nB].first,node_ends[nB].second.first-node_ends[nA].second.second+1);
        }
    }
    return linkage;
}

DistanceGraph LinkageMaker::make_longreads_multilinkage(const LongReadsMapper &lorm, uint64_t min_map_size,
                                                       float min_map_id, bool real_read_size, int32_t unmapped_end) {
    check_selected_nodes();
    auto filtered_read_mappings=lorm.improve_mappings(lorm.filter_mappings_by_size_and_id(min_map_size,min_map_id));
    DistanceGraph ldg(dg.sdg);
    std::vector<Link> linkage;
    auto lormidx=0;
    for (; lormidx < dg.sdg.ws.long_reads_datastores.size(); ++lormidx){
        if (lorm.datastore.filename == dg.sdg.ws.long_reads_datastores[lormidx].filename) break;
    }
    //for each read's filtered mappings:
    for(int64_t rid=0;rid<filtered_read_mappings.size();++rid) {
        if (filtered_read_mappings[rid].empty()) continue;
        auto newlinks=mappings_to_multilinkage(filtered_read_mappings[rid],(real_read_size ? lorm.datastore.read_to_fileRecord[rid].record_size : 0), unmapped_end);
        for (auto &l:newlinks){
            l.support.type=SupportType::LongRead;
            l.support.index=lormidx;
            l.support.id=rid;
        }
        linkage.insert(linkage.end(),newlinks.begin(),newlinks.end());
    }
    for (auto l:linkage) {
        ldg.add_link(l.source,l.dest,l.dist,l.support);
    }
    return ldg;
}

//DistanceGraph LinkageMaker::make_longreads_multilinkage(const std::string &longreads_datastore_name, uint64_t min_map_size,
//                                                        float min_map_id, bool real_read_size, int32_t unmapped_end) {
//    const auto& lorm = dg.sdg.ws.get_long_reads_datastore(longreads_datastore_name).mapper;
//    make_longreads_multilinkage(lorm, min_map_size, min_map_id, real_read_size, unmapped_end);
//}



DistanceGraph LinkageMaker::make_paired10x_multilinkage(const PairedReadsMapper &prm, const LinkedReadsMapper &lirm, float min_tnscore, bool fr,
                                                        uint64_t read_offset, int64_t min_dist_size) {
    check_selected_nodes();
    uint16_t prmidx=0;
    for (; prmidx < dg.sdg.ws.paired_reads_datastores.size(); ++prmidx){
        if (prm.datastore.filename == dg.sdg.ws.paired_reads_datastores[prmidx].filename) break;
    }
    DistanceGraph ldg(dg.sdg);
    uint64_t unmapped(0),same(0),non_neighbours(0),used(0);

    std::vector<int32_t> read_first_pos(prm.read_direction_in_node.size());
    for (auto &nm:prm.reads_in_node){
        for (auto &m:nm) read_first_pos[m.read_id]=m.first_pos;
    }
    auto sdist=prm.size_distribution();
    int64_t isize=0;
    uint64_t max_count=0;
    if (min_dist_size==-1){
        if (fr) min_dist_size=0;
        else min_dist_size=1000;
    }
    for (uint64_t i=min_dist_size;i<sdist.size();++i){
        if (sdist[i]>max_count) {
            isize=i;
            max_count=sdist[i];
        }
    }

    for (uint64_t rid1=1; rid1<prm.read_to_node.size()-1; rid1+=2){
        uint64_t rid2=rid1+1;
        sgNodeID_t n1=prm.read_to_node[rid1];
        sgNodeID_t n2=prm.read_to_node[rid2];
        if (n1==0 or n2==0) {
            ++unmapped;
            continue;
        }
        if (n1==n2) {
            ++same;
            continue;
        }
        //Check neighbourhod (any direction passes)
        bool tnpass=false;
        for (auto &tn:lirm.tag_neighbours[n1]) {
            if (tn.node == n2 and tn.score >= min_tnscore) {
                tnpass = true;
                break;
            };
        }
        if (not tnpass) {
            for (auto &tn:lirm.tag_neighbours[n2]) {
                if (tn.node == n1 and tn.score >= min_tnscore) {
                    tnpass = true;
                    break;
                };
            }
        }
        if (not tnpass) {
            ++non_neighbours;
            continue;
        }
        //orient nodes as per connection ends
        //if (fr!=prm.read_direction_in_node[rid1]) n1=-n1;
        //if (fr!=prm.read_direction_in_node[rid2]) n2=-n2;

        int d1,d2;
        if (fr==prm.read_direction_in_node[rid1]){
            d1=read_first_pos[rid1];
            n1=-n1;
        }
        else d1=dg.sdg.nodes[n1].sequence.size()-read_first_pos[rid1];
        if (fr==prm.read_direction_in_node[rid2]) {
            d2=read_first_pos[rid2];
            n2=-n2;
        }
        else d2=dg.sdg.nodes[n2].sequence.size()-read_first_pos[rid2];

        ldg.add_link(n1,n2,isize-d1-d2,{SupportType::PairedRead,prmidx,rid1});
        ++used;
    }
    sdglib::OutputLog()<<"From lmp10x: "<<unmapped<<" unmapped,  "<<same<<" same,  "<<non_neighbours<<" non-neighbours,  "<<used<<" used"<<std::endl;
    return ldg;
}

/**
 * This goes read by read, and filters the mappings by finding a set of linked nodes that maximises 1-cov of the read
 *
 * Unfiltered mappings read from mappings and results stored in filtered_read_mappings, which is cleared.
 *
 * @param lrm a LinkedReadMapper with mapped reads, over the same graph this mapper has mapped Long Reads.
 * @param min_size minimum size of the read to filter mappings.
 * @param min_tnscore minimum neighbour score on linked reads
 */
std::vector<std::vector<LongReadMapping>> LinkageMaker::filter_mappings_with_linked_reads(const LongReadsMapper &lorm, const LinkedReadsMapper &lrm, uint32_t min_size,  float min_tnscore) {
    if (lrm.tag_neighbours.empty()) {
        sdglib::OutputLog()<<"Can't filter mappings because there are no tag_neighbours on the LinkedReadsMapper"<<std::endl;
        return {};
    }
    sdglib::OutputLog()<<"Filtering mappings"<<std::endl;
    std::vector<std::vector<LongReadMapping>> filtered_read_mappings;
    filtered_read_mappings.resize(lorm.datastore.size()+1);
    std::atomic<uint64_t> rshort(0), runmapped(0), rnoresult(0), rfiltered(0);
#pragma omp parallel
    {
        LongReadHaplotypeMappingsFilter hap_filter(lorm,lrm);
        //run hap_filter on every read

#pragma omp for schedule(static, 10)
        for (uint64_t rid = 0; rid <= lorm.datastore.size(); ++rid) {
            hap_filter.set_read(rid);
            if (hap_filter.read_seq.size()<min_size){
                ++rshort;
                continue;
            }
            if (hap_filter.mappings.empty()) {
                ++runmapped;
                continue; // no mappings to group
            }
            hap_filter.generate_haplotypes_from_linkedreads(min_tnscore);
            if (hap_filter.haplotype_scores.empty()) {
                ++rnoresult;
                continue; // no haplotypes to score (all nodes too short)
            }
            hap_filter.score_coverage(1);
            hap_filter.score_window_winners(3);
            std::sort(hap_filter.haplotype_scores.rbegin(),hap_filter.haplotype_scores.rend());
            if (hap_filter.haplotype_scores[0].score>1.5) {
                std::set<sgNodeID_t> winners;
                for (auto &hn:hap_filter.haplotype_scores[0].haplotype_nodes) winners.insert(hn);
                for (auto &m:hap_filter.mappings) if (winners.count(llabs(m.node))) filtered_read_mappings[rid].emplace_back(m);
                ++rfiltered;
            } else ++rnoresult;
        }
    }
    //sdglib::OutputLog()<<"too short: "<<rshort<<"  no mappings: "<<runmapped<<"  no winner: "<<rnoresult<<"  filtered: "<<rfiltered<<std::endl;
    return filtered_read_mappings;
}

DistanceGraph LinkageMaker::make_long10x_multilinkage(const LongReadsMapper &lorm, const LinkedReadsMapper &lrm, uint32_t min_size,  float min_tnscore, bool real_read_size, int32_t unmapped_end) {
    check_selected_nodes();
    auto filtered_read_mappings=filter_mappings_with_linked_reads(lorm,lrm,min_size,min_tnscore);
    DistanceGraph ldg(dg.sdg);
    std::vector<Link> linkage;
    auto lormidx=0;
    for (; lormidx < dg.sdg.ws.long_reads_datastores.size(); ++lormidx){
        if (lorm.datastore.filename == dg.sdg.ws.long_reads_datastores[lormidx].filename) break;
    }
    //for each read's filtered mappings:
    for(int64_t rid=0;rid<filtered_read_mappings.size();++rid) {
        if (filtered_read_mappings[rid].empty()) continue;
        auto newlinks=mappings_to_multilinkage(filtered_read_mappings[rid],(real_read_size ? lorm.datastore.read_to_fileRecord[rid].record_size : 0), unmapped_end);
        for (auto &l:newlinks){
            l.support.type=SupportType::LongRead;
            l.support.index=lormidx;
            l.support.id=rid;
        }
        linkage.insert(linkage.end(),newlinks.begin(),newlinks.end());
    }
    for (auto l:linkage) {
        ldg.add_link(l.source,l.dest,l.dist,l.support);
    }
    return ldg;
}

DistanceGraph LinkageMaker::make_long10x_multilinkage(const std::string &lorm_name, const std::string &lrm_name, uint32_t min_size,  float min_tnscore, bool real_read_size, int32_t unmapped_end) {
    const auto& lorm = dg.sdg.ws.get_long_reads_datastore(lorm_name).mapper;
    const auto& lrm = dg.sdg.ws.get_linked_reads_datastore(lrm_name).mapper;

    make_long10x_multilinkage(lorm, lrm, min_size, min_tnscore, real_read_size, unmapped_end);
}

void LinkageMaker::check_selected_nodes() {
    if (std::find(selected_nodes.begin(), selected_nodes.end(), true) == selected_nodes.cend()) {
        throw std::runtime_error("No nodes selected for multilinkage");
    }
}
