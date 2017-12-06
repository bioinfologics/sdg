//
// Created by Bernardo Clavijo (EI) on 10/11/2017.
//

#include "Scaffolder.hpp"

void Scaffolder::pop_unsupported_shortbubbles() {
    std::cout<<"popping unsupported bubbles!"<<std::endl;
    uint64_t max_length=1000;
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

void Scaffolder::find_canonical_repeats(){
    const int required_support=3;
    uint64_t count=0, l700=0,l2000=0,l4000=0,l10000=0,big=0,checked=0,solvable=0;

    for (sgNodeID_t n=1;n<sg.nodes.size();++n){
        if (sg.get_fw_links(n).size()==2 and sg.get_bw_links(n).size()==2){
            ++count;
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
//
//                    std::cout << s11 << " " << s22 << " / " << s21 << " " << s12 << " " << std::endl;
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