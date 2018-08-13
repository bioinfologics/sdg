//
// Created by Bernardo Clavijo (EI) on 25/05/2018.
//

#include <sglib/logger/OutputLog.h>
#include "LinkageDiGraph.hpp"
void LinkageDiGraph::add_link(sgNodeID_t source, sgNodeID_t dest, int32_t d) {
    if (llabs(source)>=links.size()) links.resize(llabs(source)+1);
    if (llabs(dest)>=links.size()) links.resize(llabs(dest)+1);
    Link l(source,dest,d);
    for (auto el:links[(source > 0 ? source : -source)]) if (el.dest==dest and el.dist==d) return;
    links[(source > 0 ? source : -source)].emplace_back(l);
    std::swap(l.source,l.dest);
    links[(dest > 0 ? dest : -dest)].emplace_back(l);
}

void LinkageDiGraph::add_links(const LinkageDiGraph &other) {
    for (auto &lv:other.links) for (auto l:lv) add_link(l.source,l.dest,l.dist);
}

void LinkageDiGraph::remove_link(sgNodeID_t source, sgNodeID_t dest) {
    if (llabs(source)>=links.size() or llabs(dest)>=links.size()) return;
    auto & slinks = links[(source > 0 ? source : -source)];
    slinks.erase(std::remove(slinks.begin(), slinks.end(), Link(source,dest,0)), slinks.end());
    auto & dlinks = links[(dest > 0 ? dest : -dest)];
    dlinks.erase(std::remove(dlinks.begin(), dlinks.end(), Link(dest,source,0)), dlinks.end());

}

void LinkageDiGraph::disconnect_node(sgNodeID_t node) {
    for (auto fwl:get_fw_links(node)) remove_link(fwl.source,fwl.dest);
    for (auto bwl:get_bw_links(node)) remove_link(bwl.source,bwl.dest);
}

std::vector<Link> LinkageDiGraph::get_fw_links( sgNodeID_t n) const {
    std::vector<Link> r;
    if (llabs(n)>=links.size()) return r;
    for (auto &l:links[(n>0 ? n : -n)]) if (l.source==-n) r.emplace_back(l);
    return r;
}

std::vector<Link> LinkageDiGraph::get_bw_links(sgNodeID_t n) const {
    return get_fw_links (-n);
}

//returns a list of all fw nodes up to radius jumps away.
std::set<sgNodeID_t> LinkageDiGraph::fw_reached_nodes(sgNodeID_t n, int radius) const {
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

void LinkageDiGraph::remove_transitive_links(int radius) {
    std::cout<<"removing transitive connections..."<<std::endl;
    std::set<std::pair<sgNodeID_t,sgNodeID_t>> indirect;
    //TODO: check mutually transitive connections (A->B, B->C, A->C, C->B)
    for (auto fn=1;fn<links.size();++fn) {
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
                        indirect.insert(std::make_pair(-n,c));
                    }
                }

            }
//            std::cout<<"Connections for node " <<labs(n) << (n>0 ? " FW":" BW")<<":"<<std::endl;
//            for (auto l:neighbours) {
//                std::cout<<l;
//                if (indirect.count(l)>0) std::cout<<" indirect"<<std::endl;
//                else std::cout<<" DIRECT!"<<std::endl;
//            }
        }

    }
    std::cout<<indirect.size()<<" transitive connections found"<<std::endl;
    for (auto ic:indirect) remove_link(ic.first,ic.second);
    std::cout<<"removing transitive connections DONE"<<std::endl;
}

void LinkageDiGraph::report_connectivity() {
    uint64_t solved=0,solved_complex=0,solved_disconnected=0,complex=0,complex_disconected=0;
    for (auto n=1;n<sg.nodes.size();++n){
        if (sg.nodes[n].status==sgNodeDeleted) continue;
        auto fc=get_fw_links(n).size();
        auto bc=get_bw_links(n).size();
        if (fc<bc) std::swap(fc,bc);
        if (bc==0) {
            if (fc==1) ++solved_disconnected;
            else if (fc>1) ++complex_disconected;
        }
        else if (bc==1){
            if (fc==1) ++solved;
            else if (fc>1) ++solved_complex;
        }
        else ++complex;
    }
    sglib::OutputLog()<<"Connected node types:  1-1: "<<solved<<"  1-0: "<<solved_disconnected<<"  1-N: "<<solved_complex
                       <<"  N-N: "<<complex<<"  N-0: "<<complex_disconected<<std::endl;
}


std::vector<std::vector<sgNodeID_t>> LinkageDiGraph::get_all_lines(uint16_t min_nodes, uint64_t min_total_size) const {
    std::vector<std::vector<sgNodeID_t>> unitigs;
    std::vector<bool> used(sg.nodes.size(),false);

    for (auto n=1;n<sg.nodes.size();++n){
        if (used[n] or sg.nodes[n].status==sgNodeDeleted) continue;
        used[n]=true;
        std::vector<sgNodeID_t> path={n};

        //two passes: 0->fw, 1->bw, path is inverted twice, so still n is +
        for (auto pass=0; pass<2; ++pass) {
            //walk til a "non-unitig" junction
            for (auto fn = get_fw_links(path.back()); fn.size() == 1; fn = get_fw_links(path.back())) {
                if (fn[0].dest != n and fn[0].dest != -n and get_bw_links(fn[0].dest).size() == 1) {
                    path.emplace_back(fn[0].dest);
                    used[fn[0].dest > 0 ? fn[0].dest : -fn[0].dest] = true;
                } else break;
            }
            std::vector<sgNodeID_t> rpath;
            for (auto rn=path.rbegin();rn!=path.rend();++rn) rpath.emplace_back(-*rn);
            path=rpath;
        }
        if (path.size()<min_nodes) continue;
        uint64_t total_size=0;
        for (auto n:path) total_size+=sg.nodes[llabs(n)].sequence.size();
        if (total_size<min_total_size) continue;
        unitigs.push_back(path);
    }
    return unitigs;
}

std::unordered_set<sgNodeID_t> LinkageDiGraph::get_connected_nodes() const {
    std::unordered_set<sgNodeID_t> r;
    for (auto n=1;n<links.size();++n) if (!links[n].empty()) r.insert(n);
    return std::move(r);
}

void LinkageDiGraph::dump_to_text(std::string filename) {
    std::ofstream of(filename);
    for (auto &lv:links) for (auto &l:lv){
        if (llabs(l.source)<=llabs(l.dest)) of<<l.source<<" "<<l.dest<<" "<<l.dist<<std::endl;
    }
}

void LinkageDiGraph::load_from_text(std::string filename) {
    std::ifstream link_file(filename);
    Link l(0,0,0);
    while (true){
        link_file>>l.source;
        link_file>>l.dest;
        link_file>>l.dist;
        if (link_file.eof()) break;
        add_link(l.source,l.dest,l.dist);
    }
}

std::vector<std::pair<sgNodeID_t,sgNodeID_t>> LinkageDiGraph::find_bubbles(uint32_t min_size,uint32_t max_size) {
    std::vector<std::pair<sgNodeID_t,sgNodeID_t>> r;
    std::vector<bool> used(sg.nodes.size(),false);
    sgNodeID_t n1,n2;
    size_t s1,s2;
    for (n1=1;n1<sg.nodes.size();++n1){
        if (used[n1]) continue;
        //get "topologically correct" bubble: prev -> [n1 | n2] -> next
        s1=sg.nodes[n1].sequence.size();
        if (s1<min_size or s1>max_size) continue;

        auto fwl=get_fw_links(n1);
        if (fwl.size()!=1) continue;
        auto bwl=get_bw_links(n1);
        if (bwl.size()!=1) continue;
        auto next=fwl[0].dest;
        auto prev=bwl[0].dest;
        auto parln=get_bw_links(next);
        if (parln.size()!=2) continue;
        auto parlp=get_fw_links(-prev);
        if (parlp.size()!=2) continue;

        if (parlp[0].dest!=n1) n2=parlp[0].dest;
        else n2=parlp[1].dest;

        if (n2!=-parln[0].dest and n2!=-parln[1].dest) continue;
        s2=sg.nodes[llabs(n2)].sequence.size();
        if (s2<min_size or s2>max_size) continue;

        used[n1]=true;
        used[llabs(n2)]=true;
        r.emplace_back(n1,n2);

    }
    return r;
}

bool LinkageDiGraph::are_connected(sgNodeID_t n1, sgNodeID_t n2) {
    for (auto &l:links[llabs(n1)]) if (l.source==n1 and l.dest==n2) return true;
    return false;
}