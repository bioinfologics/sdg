//
// Created by Bernardo Clavijo (EI) on 10/11/2017.
//

#include "Scaffolder.hpp"

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

            if (sg.nodes[n].sequence.size()<4000) {
                // check links around the sides
                sgNodeID_t pn1=sg.get_bw_links(n)[0].dest;
                sgNodeID_t pn2=sg.get_bw_links(n)[1].dest;
                sgNodeID_t nn1=sg.get_fw_links(n)[0].dest;
                sgNodeID_t nn2=sg.get_fw_links(n)[1].dest;
                if (pn1!=pn2 and pn1!=nn1 and pn1!=nn2 and pn2!=nn1 and pn2!=nn2 and nn1!=nn2
                    and sg.nodes[pn1].sequence.size()>500
                        and sg.nodes[pn2].sequence.size()>500
                            and sg.nodes[nn1].sequence.size()>500
                                and sg.nodes[nn2].sequence.size()>500
                        ) {
                    auto s11 = count_reads_linking(pn1, nn1);
                    auto s12 = count_reads_linking(pn1, nn2);
                    auto s21 = count_reads_linking(pn2, nn1);
                    auto s22 = count_reads_linking(pn2, nn2);
                    ++checked;
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

uint64_t Scaffolder::count_reads_linking(sgNodeID_t source, sgNodeID_t dest) {
    //std::cout<<"checking support for "<<source<<" -> "<<dest<<std::endl;
    uint64_t c=0;
    for (auto &m:rmappers) {
        for (auto &rm:m.reads_in_node[(source > 0 ? source : -source)]) {
            auto prm=(rm.read_id%2 ? rm.read_id+1:rm.read_id-1);
            if (m.read_to_node[prm]==dest || m.read_to_node[prm]==-dest) ++c;//TODO:consider correct orientation and distance
            //std::cout<<rm.read_id<<"("<<rm.node<<")"<<" -> "<<prm<<"("<<m.read_to_node[prm]<<")"<<std::endl;
        }
    }
    return c;
}