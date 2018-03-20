//
// Created by Bernardo Clavijo (EI) on 10/11/2017.
//

#include "Scaffolder.hpp"

void Scaffolder::pop_unsupported_shortbubbles(uint64_t max_length) {
    std::cout<<"popping unsupported bubbles!"<<std::endl;
    //find short nodes with parallel nodes short like them
    for (sgNodeID_t n=1; n<sg.nodes.size(); ++n) {
        if (sg.nodes[n].sequence.size()>max_length) continue;
        //one link on each end, not to self, not to the same node
        if (sg.links[n].size()!=2) continue;
        if (sg.links[n][0].source==sg.links[n][1].source) continue;
        if (sg.links[n][0].dest==n or sg.links[n][0].dest==-n or sg.links[n][1].dest==n or sg.links[n][1].dest==-n) continue;
        if (sg.links[n][0].dest==sg.links[n][1].dest or sg.links[n][0].dest==-sg.links[n][1].dest) continue;
        //std::cout<<"Node "<<sg.oldnames[n]<<" passed first filter"<<std::endl;
        for (sgNodeID_t m=n+1; m<sg.nodes.size(); ++m){
            if (sg.nodes[m].sequence.size()>max_length) continue;
            //one link on each end, not to self, not to the same node
            if (sg.links[m].size()!=2) continue;
            if (sg.links[m][0].source==sg.links[m][1].source) continue;
            if (sg.links[m][0].dest==m or sg.links[m][0].dest==-m or sg.links[m][1].dest==m or sg.links[m][1].dest==-m) continue;
            if (sg.links[m][0].dest==sg.links[m][1].dest or sg.links[m][0].dest==-sg.links[m][1].dest) continue;
            if ( (sg.links[n][0].dest==sg.links[m][0].dest and sg.links[n][1].dest==sg.links[m][1].dest)
                 or (sg.links[n][1].dest==sg.links[m][0].dest and sg.links[n][0].dest==sg.links[m][1].dest)){
                //std::cout<<"found potential short bubble between #"<<n<<" "<<sg.oldnames[n]<<" and #"<<m<<" "<<sg.oldnames[m]<<std::endl;
                auto kcic1=kci.compute_compression_for_node(sg.links[n][0].dest);
                auto kcic2=kci.compute_compression_for_node(sg.links[n][1].dest);
                auto kcin=kci.compute_compression_for_node(n);
                auto kcim=kci.compute_compression_for_node(m);
                //std::cout << "kci: " << kcic1 << " -> [ " << kcin << " | " << kcim << " ] -> " << kcic1 << std::endl;
                //TODO: de-hardcode the 20% coverage rule
                if (kcic1>.8 and kcic1<1.2 and kcic2>.8 and kcic2<1.2){ //connections are unique
                    if (kcim<.2 and kcin>.5 and kcin<1.2) {
                        sg.remove_node(m);
                        std::cout<<"removed #"<<m<<std::endl;
                    }
                    if (kcin<.2 and kcim>.5 and kcim<1.2) {
                        sg.remove_node(n);
                        std::cout<<"removed #"<<n<<std::endl;
                    }

                }
            }
        }
    }
}

// TODO: Adapt for repeat resolution.
void Scaffolder::find_canonical_repeats(){
    const int required_support=3;
    uint64_t count=0, l700=0,l2000=0,l4000=0,l10000=0,big=0,checked=0,solvable=0;

    for (sgNodeID_t n=1;n<sg.nodes.size();++n){
        if (sg.get_fw_links(n).size()>=2 and sg.get_bw_links(n).size()>=2){
            ++count;

            //std::cout << "Considering node " << n << std::endl;

            if (sg.nodes[n].sequence.size()<700) ++l700;
            else if (sg.nodes[n].sequence.size()<2000) ++l2000;
            else if (sg.nodes[n].sequence.size()<4000) ++l4000;
            else if (sg.nodes[n].sequence.size()<10000) ++l10000;
            else ++big;
            // check links around the sides
            sgNodeID_t pn1=sg.get_bw_links(n)[0].dest;
            sgNodeID_t pn2=sg.get_bw_links(n)[1].dest;
            sgNodeID_t nn1=sg.get_fw_links(n)[0].dest;
            sgNodeID_t nn2=sg.get_fw_links(n)[1].dest;
            if (pn1==-pn2 or nn1==-nn2) {
                auto l=pn1;
                if (pn1==-pn2) l=(pn1>0?pn1:pn2);
                else l=(nn1>0?nn1:nn2);
                std::cout<<"Pin-hole loop found at old "<<sg.oldnames[n]<<" ("<<sg.nodes[n].sequence.size()
                         <<"bp) around edge "<<sg.oldnames[l]<<" ("<<sg.nodes[l].sequence.size() <<"bp)"<<std::endl;
            }

            if (sg.nodes[n].sequence.size()<4000) {

                if (pn1!=pn2 and pn1!=nn1 and pn1!=nn2 and pn2!=nn1 and pn2!=nn2 and nn1!=nn2
                    and sg.nodes[(pn1>0?pn1:-pn1)].sequence.size()>500
                        and sg.nodes[(pn2>0?pn2:-pn2)].sequence.size()>500
                            and sg.nodes[(nn1>0?nn1:-nn1)].sequence.size()>500
                                and sg.nodes[(nn2>0?nn2:-nn2)].sequence.size()>500
                        ) {
                    std::cout<<"evaluating trivial repeat at "<<n<<"("<<sg.nodes[n].sequence.size()<<"bp)"<<std::endl;
                    std::cout<<"PREV: "<<pn1<<"("<<sg.nodes[(pn1>0?pn1:-pn1)].sequence.size()<<"bp) "<<pn2<<"("<<sg.nodes[(pn2>0?pn2:-pn2)].sequence.size()<<"bp) "<<std::endl;
                    std::cout<<"NEXT: "<<nn1<<"("<<sg.nodes[(nn1>0?nn1:-nn1)].sequence.size()<<"bp) "<<nn2<<"("<<sg.nodes[(nn2>0?nn2:-nn2)].sequence.size()<<"bp) "<<std::endl;
                    auto s11 = count_reads_linking(pn1, nn1);
                    auto s12 = count_reads_linking(pn1, nn2);
                    auto s21 = count_reads_linking(pn2, nn1);
                    auto s22 = count_reads_linking(pn2, nn2);
                    ++checked;
//                    std::cout<<"Links for PREV#1 ("<<pn1<<")"<<std::endl;
//                    for (unsigned li=0; li<rmappers.size();++li)
//                        for (auto &rl:all_read_links(pn1,li)) std::cout<<"LIB"<<li<<": "<<rl.first<<"("<<rl.second<<")"<<std::endl;
//                    std::cout<<"Links for PREV#2 ("<<pn2<<")"<<std::endl;
//                    for (unsigned li=0; li<rmappers.size();++li)
//                        for (auto &rl:all_read_links(pn2,li)) std::cout<<"LIB"<<li<<": "<<rl.first<<"("<<rl.second<<")"<<std::endl;

                    std::cout << s11 << " " << s22 << " / " << s21 << " " << s12 << " " << std::endl;
                    //check reads supporting pn1->nn1 and pn2->nn2
                    if (s11 >= required_support and s22 >= required_support and s21 < required_support and
                        s12 < required_support)
                        ++solvable;
                    //check reads supporting pn1->nn2 and pn2->nn1
                    if (s11 < required_support and s22 < required_support and s21 >= required_support and
                        s12 >= required_support)
                        ++solvable;
                }

            }
        }
    }
    std::cout<<"Candidates for canonical repeat expansion:                    "<<count<<std::endl;
    std::cout<<"Candidates for canonical repeat expansion <700bp:             "<<l700<<std::endl;
    std::cout<<"Candidates for canonical repeat expansion >700bp & <2000bp:   "<<l2000<<std::endl;
    std::cout<<"Candidates for canonical repeat expansion >2000bp & <4000bp:  "<<l4000<<std::endl;
    std::cout<<"Candidates for canonical repeat expansion >4000bp & <10000bp: "<<l10000<<std::endl;
    std::cout<<"Candidates for canonical repeat expansion >10000bp:           "<<big<<std::endl;
    std::cout<<"Trivially solvable canonical repeats:                         "<<solvable<<"/"<<checked<<std::endl;
}

std::vector<SequenceSubGraph> Scaffolder::get_all_bubbly_subgraphs(uint32_t maxsubgraphs) {
    std::vector<SequenceSubGraph> subgraphs;
    std::vector<bool> used(sg.nodes.size(),false);
    const double min_c1=0.75,max_c1=1.25,min_c2=1.5,max_c2=2.5;
    /*
     * the loop always keep the first and the last elements as c=2 collapsed nodes.
     * it starts with a c=2 node, and goes thorugh all bubbles fw, then reverts the subgraph and repeats
     */
    SequenceSubGraph subgraph(sg);
    for (auto n=1;n<sg.nodes.size();++n){
        if (used[n] or sg.nodes[n].status==sgNodeDeleted) continue;
        auto frontkci=kci.compute_compression_for_node(n);
        if (frontkci>max_c2 or frontkci<min_c2) continue;
        used[n]=true;
        subgraph.nodes.clear();

        subgraph.nodes.push_back(n);

        //two passes: 0->fw, 1->bw, path is inverted twice, so still n is +
        for (auto pass=0; pass<2; ++pass) {
            //while there's a possible bubble fw.
            for (auto fn = sg.get_fw_links(subgraph.nodes.back()); fn.size() == 2; fn = sg.get_fw_links(subgraph.nodes.back())) {
                //if it is not a real bubble, get out.
                if (sg.get_bw_links(fn[0].dest).size()!=1 or sg.get_bw_links(fn[0].dest).size()!=1) break;
                auto fl1=sg.get_fw_links(fn[0].dest);
                if (fl1.size()!=1) break;
                auto fl2=sg.get_fw_links(fn[1].dest);
                if (fl2.size()!=1) break;
                if (fl2[0].dest!=fl1[0].dest) break;

                auto p1kci=kci.compute_compression_for_node(fn[0].dest);
                if (p1kci<min_c1 or p1kci>max_c1) break;
                auto p2kci=kci.compute_compression_for_node(fn[1].dest);
                if (p2kci<min_c1 or p2kci>max_c1) break;

                auto next_end=fl2[0].dest;
                auto next_end_kci=kci.compute_compression_for_node(next_end);
                if (next_end_kci<min_c2 or next_end_kci>max_c2) break;

                //all conditions met, update subgraph
                subgraph.nodes.push_back(fn[0].dest);
                subgraph.nodes.push_back(fn[1].dest);
                subgraph.nodes.push_back(next_end);
                used[(next_end>0?next_end:-next_end)]=true;

            }
            SequenceSubGraph new_subgraph(sg);
            for (auto it=subgraph.nodes.rbegin();it<subgraph.nodes.rend();++it) new_subgraph.nodes.push_back(-*it);
            std::swap(new_subgraph.nodes,subgraph.nodes);
        }
        if (subgraph.nodes.size()>6) {
            subgraphs.push_back(subgraph);
            if (subgraphs.size()==maxsubgraphs) break;
            //std::cout<<"Bubbly path found: ";
            //for (auto &n:subgraph) std::cout<<"  "<<n<<" ("<<sg.nodes[(n>0?n:-n)].sequence.size()<<"bp)";
            //std::cout<<std::endl;
        }
    }
    return subgraphs;
}

void Scaffolder::expand_bubbly_subgraphs() {
    std::cout<<"Finding bubbly subgraphs"<<std::endl;
    auto bubbly_subgraphs=get_all_bubbly_subgraphs();
    std::cout<<"Expanding bubbly subgraphs"<<std::endl;
    for (auto i=1;i<rmappers[0].read_to_tag.size();++i){
        if(rmappers[0].read_to_tag[i]==0) {
            if (i%2==1) rmappers[0].read_to_tag[i]=rmappers[0].read_to_tag[i+1];
            else rmappers[0].read_to_tag[i]=rmappers[0].read_to_tag[i-1];
        }
    }
    for (auto bubblysg:bubbly_subgraphs){
        if (bubblysg.nodes.size()<16) continue; //skip small bubbly paths
        //get all tags involved in region projected into nodes in the region
        std::cout<<std::endl<<"Analysing bubbly subgraph:"<<std::endl;

        std::set<uint32_t> tags;
        for (auto node:bubblysg.nodes){
            auto n=(node>0?node:-node);
            for (auto r:rmappers[0].reads_in_node[n]) tags.insert(rmappers[0].read_to_tag[r.read_id]);
        }
        std::cout<<"There are "<<tags.size()<<" tags in the region"<<std::endl;
        //create the node/tag read matrix.
        std::map<uint32_t,uint32_t> tags_to_local;
        std::vector<uint32_t>local_to_tag;
        {
            auto i=0;
            for (auto it = tags.begin(); it != tags.end(); ++it, ++i) {
                tags_to_local[*it] = i;
                local_to_tag.push_back(*it);
            }
        }
        std::cout<<"Local tag direct and reverse indexes done"<<std::endl;
        std::vector<std::vector<uint64_t>> node_tag_readcount;
        for (auto node:bubblysg.nodes) {
            node_tag_readcount.emplace_back(std::vector<uint64_t>(tags.size()));
            auto n = (node > 0 ? node : -node);
            std::cout << "processing node #" << n << ":" << node;
            for (auto r:rmappers[0].reads_in_node[n]) ++node_tag_readcount.back()[tags_to_local[rmappers[0].read_to_tag[r.read_id]]];
        }
        std::cout << "dumping matrix"<<std::endl;
        std::ofstream node_tag_mx_file("node_tag_mx_"+std::to_string(bubblysg.nodes[0])+".csv");
        node_tag_mx_file<<"tag";
        for (auto n:bubblysg.nodes)node_tag_mx_file<<",node_"<<n;
        node_tag_mx_file<<std::endl;
        for (auto i=0;i<tags.size();++i){
            node_tag_mx_file<<local_to_tag[i];
            for (auto j=0;j<bubblysg.nodes.size();++j) node_tag_mx_file<<","<<node_tag_readcount[j][i];
            node_tag_mx_file<<std::endl;
        }


        //find the right node groups and create two paths


        //expand the bubles using sg.join_path

    }
}

std::vector<std::pair<sgNodeID_t,sgNodeID_t>> Scaffolder::get_all_haplotype_pairs(uint32_t maxpairs) {
    std::vector<std::pair<sgNodeID_t,sgNodeID_t>> hps;
    std::vector<bool> used(sg.nodes.size(),false);
    //TODO: check the coverages are actually correct?
    const double min_c1=0.7,max_c1=1.30,min_c2=1.8,max_c2=2.8;
    /*
     * the loop always keep the first and the last elements as c=2 collapsed nodes.
     * it starts with a c=2 node, and goes thorugh all bubbles fw, then reverts the subgraph and repeats
     */
    std::pair<sgNodeID_t,sgNodeID_t> hap={0,0};
    for (auto n=1;n<sg.nodes.size();++n){
        if (sg.nodes[n].status==sgNodeDeleted) continue;
        auto frontkci=kci.compute_compression_for_node(n);
        if (frontkci>max_c2 or frontkci<min_c2) continue;
        auto m=n;
        //two passes: 0->fw, 1->bw,
        for (auto pass=0; pass<2; ++pass,m=-m) {

            //check bubble going forward ------------
            auto fw_l = sg.get_fw_links(m);
            //fork opening
            if (fw_l.size() != 2) continue;
            hap.first = fw_l[0].dest;
            hap.second = fw_l[1].dest;
            if (hap.first == n or hap.first == -n or hap.second == n or hap.second == -n or hap.first == hap.second) continue;
            //fork.closing
            auto hap0f = sg.get_fw_links(hap.first);
            auto hap1f = sg.get_fw_links(hap.second);
            if (hap0f.size() != 1 or hap1f.size() != 1 or hap0f[0].dest != hap1f[0].dest) continue;
            auto h0kc = kci.compute_compression_for_node(hap.first);
            if (h0kc > max_c1 or h0kc < min_c1) continue;
            auto h1kc = kci.compute_compression_for_node(hap.second);
            if (h1kc > max_c1 or h1kc < min_c1) continue;
            auto ekc = kci.compute_compression_for_node(hap0f[0].dest);
            if (ekc > max_c2 or ekc < min_c2) continue;
            if (used[(hap.first>0?hap.first:-hap.first)] or used[(hap.second>0?hap.second:-hap.second)]){
                if (not used[(hap.second>0?hap.second:-hap.second)] or not used[(hap.second>0?hap.second:-hap.second)])
                    std::cout<<"WARNING: unusual use pattern on HSPNP"<<std::endl;
                continue;
            }
            hps.push_back(hap);
            if (hps.size()%100==0) std::cout<<hps.size()<<" haplotype pairs found"<<std::endl;
            used[(hap.first>0?hap.first:-hap.first)]=true;
            used[(hap.second>0?hap.second:-hap.second)]=true;
        }
    }
    return hps;
}

uint64_t Scaffolder::count_reads_linking(sgNodeID_t source, sgNodeID_t dest) {
    //std::cout<<"checking support for "<<source<<" -> "<<dest<<std::endl;
    uint64_t c=0;
    for (auto &m:rmappers) {
        for (auto &rm:m.reads_in_node[(source > 0 ? source : -source)]) {
            auto prm=(rm.read_id%2 ? rm.read_id+1:rm.read_id-1);
            if (m.read_to_node[prm]==dest || m.read_to_node[prm]==-dest) ++c;// \todo consider correct orientation and distance
            //std::cout<<rm.read_id<<"("<<rm.node<<")"<<" -> "<<prm<<"("<<m.read_to_node[prm]<<")"<<std::endl;
        }
    }
    return c;
}

std::vector<std::pair<sgNodeID_t,uint64_t>> Scaffolder::all_read_links(sgNodeID_t source, unsigned lib) {
    std::vector<std::pair<sgNodeID_t,uint64_t>> links;


    for (auto &rm:rmappers[lib].reads_in_node[(source > 0 ? source : -source)]) {
        auto prm=(rm.read_id%2 ? rm.read_id+1:rm.read_id-1);
        auto dest=rmappers[lib].read_to_node[prm];
        if (dest==source or dest==-source) continue;
        if (!rm.rev) dest=-dest;//TODO:check this, configure library direction

        auto it = links.begin();
        for (; it != links.end(); ++it) if(it->first==dest) break;
        if(it == links.end()) links.emplace_back(dest,1);
        else ++(it->second);
    }
    std::sort(links.begin(),links.end(),[](std::pair<sgNodeID_t,uint64_t> &left, std::pair<sgNodeID_t,uint64_t> &right) {
        return left.second > right.second;
    });
    return links;
}