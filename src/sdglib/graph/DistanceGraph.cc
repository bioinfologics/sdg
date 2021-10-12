//
// Created by Bernardo Clavijo (EI) on 25/05/2018.
//

#include "DistanceGraph.hpp"
#include <fstream>
#include <cmath>
#include <sdglib/utilities/OutputLog.hpp>
#include "SequenceDistanceGraph.hpp"
#include <sdglib/views/NodeView.hpp>
#include <sdglib/views/TangleView.hpp>
#include <sstream>
#include <sdglib/workspace/WorkSpace.hpp>
#include <numeric>

void DistanceGraph::add_link(sgNodeID_t source, sgNodeID_t dest, int32_t d, Support support) {
    if (llabs(source)>=links.size()) links.resize(llabs(source)+1);
    if (llabs(dest)>=links.size()) links.resize(llabs(dest)+1);
    Link l(source,dest,d,support);
    links[llabs(source)].emplace_back(l);
    if (l.source!=l.dest) {
        std::swap(l.source, l.dest);
        links[llabs(dest)].emplace_back(l);
    }
}

void DistanceGraph::copy_links(const DistanceGraph &other) {
    for (auto &lv:other.links) for (auto l:lv) add_link(l.source,l.dest,l.dist,l.support);
}

bool DistanceGraph::remove_link(sgNodeID_t source, sgNodeID_t dest) {
    auto & slinks = links[(source > 0 ? source : -source)];
    auto slinksLen = slinks.size();
    slinks.erase(std::remove_if(slinks.begin(), slinks.end(), [source,dest] (auto l) {return l.source==source and l.dest==dest;}), slinks.end());
    auto & dlinks = links[(dest > 0 ? dest : -dest)];
    auto dlinksLen = dlinks.size();
    dlinks.erase(std::remove_if(dlinks.begin(), dlinks.end(), [source,dest] (auto l) {return l.source==dest and l.dest==source;}), dlinks.end());
    if (slinks.size() != slinksLen or dlinks.size() != dlinksLen) return true;
    return false;
}

bool DistanceGraph::remove_link(sgNodeID_t source, sgNodeID_t dest, int32_t d, Support s) {
    auto & slinks = links[(source > 0 ? source : -source)];
    auto slinksLen = slinks.size();
    slinks.erase(std::remove_if(slinks.begin(), slinks.end(), [source,dest,d,s](const Link &l) {
        return std::tie(source,dest,d,s)==std::tie(l.source,l.dest,l.dist,l.support);
    }), slinks.end());
    auto & dlinks = links[(dest > 0 ? dest : -dest)];
    auto dlinksLen = dlinks.size();
    dlinks.erase(std::remove_if(dlinks.begin(), dlinks.end(), [source,dest,d,s](const Link &l) {
        return std::tie(source,dest,d,s)==std::tie(l.dest,l.source,l.dist,l.support);
    }), dlinks.end());
    if (slinks.size() != slinksLen or dlinks.size() != dlinksLen) return true;
    return false;
}

void DistanceGraph::disconnect_node(sgNodeID_t node) {
    for (auto fwl:get_fw_links(node)) remove_link(fwl.source,fwl.dest);
    for (auto bwl:get_bw_links(node)) remove_link(bwl.source,bwl.dest);
}

std::vector<Link> DistanceGraph::get_fw_links( sgNodeID_t n) const {
    std::vector<Link> r;
    if (llabs(n)>=links.size()) return r;
    for (auto &l:links[(n>0 ? n : -n)]) if (l.source==-n) r.emplace_back(l);
    return r;
}

std::vector<Link> DistanceGraph::get_bw_links(sgNodeID_t n) const {
    return get_fw_links (-n);
}

Link DistanceGraph::get_link(sgNodeID_t source, sgNodeID_t dest) {
    for (auto l:links[(source > 0 ? source : -source)]) {
        if (l.source==source,l.dest==dest) return l;
    }
    return Link(0,0,0);
}

std::vector<sgNodeID_t> DistanceGraph::get_next_nodes(sgNodeID_t n) const {
    std::vector<sgNodeID_t > r;
    for (auto&  l:links[llabs(n)]) if (l.source==-n) r.emplace_back(l.dest);
    return r;
}

std::vector<sgNodeID_t> DistanceGraph::get_prev_nodes(sgNodeID_t n) const {
    std::vector<sgNodeID_t > r;
    for (auto&  l:links[llabs(n)]) if (l.source==n) r.emplace_back(-l.dest);
    return r;
}

//returns a list of all fw nodes up to radius jumps away.
std::set<sgNodeID_t> DistanceGraph::fw_reached_nodes(sgNodeID_t n, int radius) const {
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

std::vector<std::pair<int,sgNodeID_t>> DistanceGraph::fw_neighbours_by_distance(sgNodeID_t n, int min_links) const {
    if (min_links<1 or llabs(n)>=links.size()) return {};
    std::map<sgNodeID_t,std::vector<int>> fndists;//neighbour and all distances to it
    std::vector<std::pair<int,sgNodeID_t>> fnd; //median distance and neighbours if with min_links
    //ndists[dest] -> [dist1, dist2, dist3...]
    for (auto &l:links[llabs(n)]){
        if (l.source==-n){
            fndists[l.dest].emplace_back(l.dist);
        }
    }
    //fnd <- [(median2, n1),(median2,n2)...]
    for (auto &fn:fndists){
        if (fn.second.size()>=min_links){
            std::sort(fn.second.begin(),fn.second.end());
            if (fn.second.size()%2==0) fnd.emplace_back(std::make_pair((fn.second[fn.second.size()/2-1]+fn.second[fn.second.size()/2])/2,fn.first));
            else fnd.emplace_back(fn.second[fn.second.size()/2],fn.first);
        }
    }
    std::sort(fnd.begin(),fnd.end());
    return fnd;
}

void DistanceGraph::remove_transitive_links(int radius) {
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

void DistanceGraph::report_connectivity() {
    uint64_t solved=0,solved_complex=0,solved_disconnected=0,complex=0,complex_disconected=0;
    for (auto n=1;n<sdg.nodes.size();++n){
        if (sdg.nodes[n].status==NodeStatus::Deleted) continue;
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
    sdglib::OutputLog()<<"Connected node types:  1-1: "<<solved<<"  1-0: "<<solved_disconnected<<"  1-N: "<<solved_complex
                       <<"  N-N: "<<complex<<"  N-0: "<<complex_disconected<<std::endl;
}


std::vector<std::vector<sgNodeID_t>> DistanceGraph::get_all_lines(uint16_t min_nodes, uint64_t min_total_size) const {
    std::vector<std::vector<sgNodeID_t>> unitigs;
    std::vector<bool> used(sdg.nodes.size(),false);

    for (auto n=1;n<sdg.nodes.size();++n){
        if (used[n] or sdg.nodes[n].status==NodeStatus::Deleted) continue;
        used[n]=true;
        std::vector<sgNodeID_t> path={n};

        //two passes: 0->fw, 1->bw, path is inverted twice, so still n is +
        for (auto pass=0; pass<2; ++pass) {
            //walk til a "non-unitig" junction
            for (auto fn = get_fw_links(path.back()); fn.size() == 1; fn = get_fw_links(path.back())) {
                if (used[llabs(fn[0].dest)]!=true and get_bw_links(fn[0].dest).size() == 1) {
                    path.emplace_back(fn[0].dest);
                    used[llabs(fn[0].dest)] = true;
                } else break;
            }
            std::vector<sgNodeID_t> rpath;
            for (auto rn=path.rbegin();rn!=path.rend();++rn) rpath.emplace_back(-*rn);
            path=rpath;
        }
        if (path.size()<min_nodes) continue;
        uint64_t total_size=0;
        for (auto n:path) total_size+=sdg.nodes[llabs(n)].sequence.size();
        if (total_size<min_total_size) continue;
        unitigs.push_back(path);
    }
    return unitigs;
}

std::unordered_set<sgNodeID_t> DistanceGraph::get_connected_nodes() const {
    std::unordered_set<sgNodeID_t> r;
    for (auto n=1;n<links.size();++n) if (!links[n].empty()) r.insert(n);
    return std::move(r);
}

std::vector<SequenceDistanceGraphPath> DistanceGraph::find_all_paths_between(sgNodeID_t from,sgNodeID_t to, int64_t max_size, int max_nodes, bool abort_on_loops) const {
    typedef struct T {
        int64_t prev;
        sgNodeID_t node;
        int node_count;
        int64_t partial_size;
        T(uint64_t a, sgNodeID_t b, int c, int64_t d) : prev(a),node(b),node_count(c),partial_size(d){};
    } pathNodeEntry_t;

    std::vector<pathNodeEntry_t> node_entries;
    node_entries.reserve(100000);
    std::vector<SequenceDistanceGraphPath> final_paths;
    std::vector<sgNodeID_t> pp;

    for(auto &fl:get_fw_links(from)) {
        if (fl.dest==to and final_paths.empty()) {
            final_paths.emplace_back(sdg,pp);
        }
        else if (sdg.nodes[llabs(fl.dest)].sequence.size()<=max_size) {
            node_entries.emplace_back(-1, fl.dest, 1, sdg.nodes[llabs(fl.dest)].sequence.size());
        }
    }


    // XXX: This loop is growing forever on node 2083801 for backbone 58, the limits are heuristics to cap the complexity of regions captured
    for (uint64_t current_index=0;current_index<node_entries.size() and node_entries.size() < 1000001 and final_paths.size() < 1001;++current_index){
        auto current_entry=node_entries[current_index];
        auto fwl=get_fw_links(current_entry.node);
        for(auto &fl:fwl) {
            if (fl.dest==to){
                //reconstruct path backwards
                pp.resize(current_entry.node_count);

                for (uint64_t ei=current_index;ei!=-1;ei=node_entries[ei].prev){
                    pp[node_entries[ei].node_count-1]=node_entries[ei].node;
                }
                //check if there are any loops
                if (abort_on_loops) {
                    for (auto i1 = 0; i1 < pp.size(); ++i1) {
                        for (auto i2 = 0; i2 < pp.size(); ++i2) {
                            if (pp[i1]==pp[i2]) {
                                return {};
                            }
                        }
                    }
                }
                //add to solutions
                final_paths.emplace_back(sdg,pp);
            }
            else {
                uint64_t new_size=current_entry.partial_size+fl.dist+sdg.nodes[llabs(fl.dest)].sequence.size();
                if (new_size<=max_size and current_entry.node_count<=max_nodes) {
                    node_entries.emplace_back(current_index, fl.dest, current_entry.node_count+1, new_size);
                }
            }
        }

    }

    bool early_exit=false;
    if (node_entries.size() == 1000000) {
        std::cout << "From " << from << " to " << to << " with max_size " << max_size << " and max_nodes " << max_nodes << " there were too many nodes!"<< std::endl;
        early_exit = true;
    }

    if (final_paths.size() == 1000) {
        std::cout << "From " << from << " to " << to << " with max_size " << max_size << " and max_nodes " << max_nodes << " there were too many paths!"<< std::endl;
        early_exit = true;
    }

    if (early_exit) return {};
    return final_paths;
}

void DistanceGraph::dump_to_text(std::string filename) {
    std::ofstream of(filename);
    for (auto &lv:links) for (auto &l:lv){
        if (llabs(l.source)<=llabs(l.dest)) of<<l.source<<" "<<l.dest<<" "<<l.dist<<std::endl;
    }
}

void DistanceGraph::load_from_text(std::string filename) {
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

std::vector<std::pair<sgNodeID_t,sgNodeID_t>> DistanceGraph::find_bubbles(uint32_t min_size,uint32_t max_size) const {
    std::vector<std::pair<sgNodeID_t,sgNodeID_t>> r;
    std::vector<bool> used(sdg.nodes.size(),false);
    sgNodeID_t n1,n2;
    size_t s1,s2;
    for (n1=1;n1<sdg.nodes.size();++n1){
        if (used[n1]) continue;
        //get "topologically correct" bubble: prev -> [n1 | n2] -> next
        s1=sdg.nodes[n1].sequence.size();
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
        s2=sdg.nodes[llabs(n2)].sequence.size();
        if (s2<min_size or s2>max_size) continue;

        used[n1]=true;
        used[llabs(n2)]=true;
        r.emplace_back(n1,n2);

    }
    return r;
}

std::vector<sgNodeID_t> DistanceGraph::find_tips(uint32_t min_size, uint32_t max_size) const {
    std::vector<sgNodeID_t> r;
    for (auto &l:links){
        if (l.size()==1 and sdg.nodes[llabs(l[0].source)].sequence.size()>=min_size and sdg.nodes[llabs(l[0].source)].sequence.size()<=max_size){
            for (auto &ol:links[llabs(l[0].dest)]){
                if (ol.source==l[0].dest and ol.dest!=l[0].source){
                    r.emplace_back(l[0].source);
                    break;
                }
            }
        }
    }
    return r;
}

std::vector<sgNodeID_t> DistanceGraph::find_self_loops(uint32_t min_size, uint32_t max_size, bool include_circles) const {
    std::vector<sgNodeID_t> r;
    for (auto &lv:links) {
        bool loop=false;
        bool other=include_circles;
        for (auto &l:lv){
            if (llabs(l.source)==llabs(l.dest)) loop=true;
            else other=true;
            if (loop and other) {
                r.emplace_back(llabs(l.source));
                break;
            }
        }
    }
    return r;
}

bool DistanceGraph::are_connected(sgNodeID_t n1, sgNodeID_t n2) const {
    if (llabs(n1)>=links.size()) return false;
    for (auto &l:links[llabs(n1)]) if (l.source==n1 and l.dest==n2) return true;
    return false;
}

uint32_t DistanceGraph::link_count(sgNodeID_t n1, sgNodeID_t n2) const {
    uint32_t c=0;
    for (auto &l:links[llabs(n1)]) if (l.source==n1 and l.dest==n2) ++c;
    return c;
}

int32_t DistanceGraph::min_distance(const sgNodeID_t & n1, const sgNodeID_t & n2) const {
    if (n1>0) {
        if (n1 >= links.size()) return INT32_MAX;
        int32_t d = INT32_MAX;
        for (const Link * l=links[n1].data();l<links[n1].data()+links[n1].size();++l)
            if (l->dest == n2 and l->source == -n1 and l->dist < d)
                d = l->dist;
        return d;
    }
    else {
        if (-n1 >= links.size()) return INT32_MAX;
        int32_t d = INT32_MAX;
        for (const Link * l=links[-n1].data();l<links[-n1].data()+links[-n1].size();++l)
            if (l->dest == n2 and l->source == -n1 and l->dist < d)
                d = l->dist;
        return d;
    }
}

void DistanceGraph::write_to_gfa1(std::string filename, const std::vector<sgNodeID_t> &selected_nodes, const std::vector<double> &depths) {
    //TODO: change gaps to sequences named gapXX without sequence, but with length
    std::unordered_set<sgNodeID_t > output_nodes(selected_nodes.begin(), selected_nodes.end());
    std::string fasta_filename;
    //check the filename ends in .gfa
    if (filename.size()>4 and filename.substr(filename.size()-4,4)==".gfa"){
        fasta_filename=filename.substr(0,filename.size()-4)+".fasta";
    }
    else throw std::invalid_argument("filename of the gfa input does not end in gfa, it ends in '"+filename.substr(filename.size()-4,4)+"'");

    std::ofstream gfaf(filename);
    if (!gfaf) throw std::invalid_argument("Can't write to gfa file");
    gfaf<<"H\tVN:Z:1.0"<<std::endl;

    std::ofstream fastaf(fasta_filename);
    if (!fastaf) throw std::invalid_argument("Can't write to fasta file");


    //load all sequences from fasta file if they're not canonical, flip and remember they're flipped
    //std::cout<<"Writing sequences to "<<fasta_filename<<std::endl;

    for (sgNodeID_t i=1;i<sdg.nodes.size();++i){
        if (sdg.nodes[i].status==NodeStatus::Deleted) continue;
        if (!output_nodes.empty() and output_nodes.count(i)==0 and output_nodes.count(-i)==0) continue;
        fastaf<<">seq"<<i<<std::endl<<sdg.nodes[i].sequence<<std::endl;
        gfaf<<"S\tseq"<<i<<"\t*\tLN:i:"<<sdg.nodes[i].sequence.size()<<"\tUR:Z:"<<fasta_filename
            <<(depths.empty() or std::isnan(depths[i])?"":"\tDP:f:"+std::to_string(depths[i]))<<std::endl;
    }

    // Go through all links and create a gap_sequence for every positive distance link
    for (const auto &ls: links) {
        for (const auto &l:ls) {
            if (l.source <= l.dest and l.dist > 0 and (output_nodes.empty() or
                                                       output_nodes.count(l.source) > 0 or output_nodes.count(-l.source) > 0 or
                                                       output_nodes.count(l.dest) > 0 or output_nodes.count(-l.dest) > 0)) {
                gfaf << "S\tgap_" << l.source << "_" << l.dest << "\t*\tLN:i:" << l.dist << std::endl;
            }
        }
    }

    for (auto &ls:links) {
        for (auto &l:ls) {
            if (l.source <= l.dest and (output_nodes.empty() or
                                        output_nodes.count(l.source) > 0 or output_nodes.count(-l.source) > 0 or
                                        output_nodes.count(l.dest) > 0 or output_nodes.count(-l.dest) > 0)) {
                if (l.dist <= 0) {
                    gfaf << "L\t";
                    if (l.source > 0) gfaf << "seq" << l.source << "\t-\t";
                    else gfaf << "seq" << -l.source << "\t+\t";

                    if (l.dest > 0) gfaf << "seq" << l.dest << "\t+\t";
                    else gfaf << "seq" << -l.dest << "\t-\t";

                    gfaf << (l.dist < 0 ? -l.dist : 0) << "M" << std::endl;
                } else {
                    gfaf << "L\t";
                    if (l.source > 0) gfaf << "seq" << l.source << "\t-\t" << "gap_" << l.source << "_" << l.dest << "\t+\t";
                    else gfaf << "seq" << -l.source << "\t+\t" << "gap_" << l.source << "_" << l.dest << "\t+\t";
                    gfaf << "0M" << std::endl;

                    gfaf << "L\t";
                    if (l.dest > 0) gfaf << "gap_" << l.source << "_" << l.dest << "\t+\t" << "seq" << l.dest << "\t+\t";
                    else gfaf << "gap_" << l.source << "_" << l.dest << "\t+\t" << "seq" << -l.dest << "\t-\t";
                    gfaf << "0M" << std::endl;
                }
            }
        }
    }
}

void DistanceGraph::write_to_gfa2(std::string filename, const std::vector<sgNodeID_t> &selected_nodes, const std::vector<double> &depths) {
    std::unordered_set<sgNodeID_t > output_nodes(selected_nodes.begin(), selected_nodes.end());
    std::string fasta_filename;
    //check the filename ends in .gfa
    if (filename.size()>4 and filename.substr(filename.size()-4,4)==".gfa"){
        fasta_filename=filename.substr(0,filename.size()-4)+".fasta";
    }
    else throw std::invalid_argument("filename of the gfa input does not end in gfa, it ends in '"+filename.substr(filename.size()-4,4)+"'");

    std::ofstream gfaf(filename);
    if (!gfaf) throw std::invalid_argument("Can't write to gfa file");
    gfaf<<"H\tVN:Z:2.0"<<std::endl;

    std::ofstream fastaf(fasta_filename);
    if (!fastaf) throw std::invalid_argument("Can't write to fasta file");


    //load all sequences from fasta file if they're not canonical, flip and remember they're flipped
    //std::cout<<"Writing sequences to "<<fasta_filename<<std::endl;

    for (sgNodeID_t i=1;i<sdg.nodes.size();++i){
        if (sdg.nodes[i].status==NodeStatus::Deleted) continue;
        if (!output_nodes.empty() and output_nodes.count(i)==0 and output_nodes.count(-i)==0) continue;
        fastaf<<">seq"<<i<<std::endl<<sdg.nodes[i].sequence<<std::endl;
        gfaf<<"S\tseq"<<i<<"\t*\tLN:i:"<<sdg.nodes[i].sequence.size()<<"\tUR:Z:"<<fasta_filename
            <<(depths.empty() or std::isnan(depths[i])?"":"\tDP:f:"+std::to_string(depths[i]))<<std::endl;
    }

    uint edge_counter(1);
    uint gap_counter(1);
    for (auto &ls:links) {
        for (auto &l:ls)
            if (l.source <= l.dest and (output_nodes.empty() or
                                        output_nodes.count(l.source) > 0 or output_nodes.count(-l.source) > 0 or
                                        output_nodes.count(l.dest) > 0 or output_nodes.count(-l.dest) > 0)) {
                if (l.dist < 0) {
                    gfaf << "E\tedge_" << edge_counter++ << "\t";
                    const auto source_size(sdg.nodes[std::abs(l.source)].sequence.size());
                    const auto dest_size(sdg.nodes[std::abs(l.dest)].sequence.size());
                    if (l.source > 0) gfaf << "seq" << l.source << "-\t";
                    else gfaf << "seq" << -l.source << "+\t";

                    if (l.dest > 0) gfaf << "seq" << l.dest << "+\t";
                    else gfaf << "seq" << -l.dest << "-\t";

                    if (l.source > 0) gfaf << source_size + l.dist << "\t" << source_size << "$\t";
                    else gfaf << "0\t" << l.dist << "\t";

                    if (l.dest > 0) gfaf << "0\t" << l.dist << "\t";
                    else gfaf << dest_size + l.dist << "\t" << dest_size << "$\t";

                    gfaf << (l.dist < 0 ? -l.dist : 0) << "M" << std::endl;
                } else {
                    gfaf << "G\tgap_"<< gap_counter++ <<"\t";
                    if (l.source > 0) gfaf << "seq" << l.source << "-\t";
                    else gfaf << "seq" << -l.source << "+\t";
                    if (l.dest > 0) gfaf << "seq" << l.dest << "+\t";
                    else gfaf << "seq" << -l.dest << "-\t";

                    gfaf << l.dist << "\t*" << std::endl;
                }
            }
    }
}

void DistanceGraph::read(std::ifstream &input_file) {
    uint64_t s;
    input_file.read((char *) &s, sizeof(s));
    name.resize(s);
    input_file.read((char *) name.data(), name.size());
    sdglib::read_flat_vectorvector(input_file, links);
}

void DistanceGraph::write(std::ofstream &output_file) {
    uint64_t s=name.size();
    output_file.write((char *) &s,sizeof(s));
    output_file.write((char *)name.data(),name.size());
    sdglib::write_flat_vectorvector(output_file, links);
}

DistanceGraph::DistanceGraph(SequenceDistanceGraph &_sdg, const std::string &_name): sdg(_sdg),name(_name) {
    if (name!="SDG") links.resize(sdg.nodes.size()); //This is a hack to allow paths in a straight-constructed DG.
}

std::string DistanceGraph::ls(int level, bool recursive) const {
    std::stringstream ss;
    std::string spacer(2*level,' ');
    uint64_t linkcount=0;
    for (auto lv:links)
        for (auto l:lv)
            if (llabs(l.source) <= llabs(l.dest)) ++linkcount;
    ss<<spacer<<"DistanceGraph "<<name<<" ("<<sdg.name<<"): "<<linkcount<<" links"<<std::endl;
    return ss.str();
}

NodeView DistanceGraph::get_nodeview(sgNodeID_t n) const {
    if (n==0 or llabs(n)>=sdg.nodes.size()) throw std::runtime_error("Trying to get a nodeview with invalid node id");
    if (sdg.nodes[llabs(n)].status==NodeStatus::Deleted) throw std::runtime_error("Trying to get a nodeview of a deleted Node");
    return NodeView(this,n);
}

std::vector<NodeView> DistanceGraph::get_all_nodeviews(bool both_directions, bool include_disconnected) const {
    uint64_t c=0;
    for (auto nidx=0;nidx<sdg.nodes.size();++nidx) {
        auto &n=sdg.nodes[nidx];
        if (n.status!=NodeStatus::Deleted and (include_disconnected or !links[nidx].empty())) {
            ++c;
            if (both_directions) ++c;
        }
    }
    std::vector<NodeView> nvs;
    nvs.reserve(c);
    for (auto nidx=0;nidx<sdg.nodes.size();++nidx) {
        auto &n=sdg.nodes[nidx];
        if (n.status!=NodeStatus::Deleted and (include_disconnected or !links[nidx].empty())) {
            nvs.emplace_back(this,nidx);
            if (both_directions) nvs.emplace_back(this,-nidx);
        }

    }
    return nvs;
}

std::vector<TangleView> DistanceGraph::get_all_anchors_tangleviews(const std::vector<sgNodeID_t> given_frontiers, bool include_disconnected) const {
//    std::cout<<"Tangles of size < "<<fsize<<" with "<<fmin_kci<<" < KCI < "<<fmax_kci<<" disconnected: "<<include_disconnected<<std::endl;
    // Make sure all elements in given frontiers are positive

    for (const auto &gf: given_frontiers){
        if (gf<=0)
            std::domain_error("All frontiers must be in the fw position (+)");
    }
    std::set<sgNodeID_t> fset(given_frontiers.begin(),given_frontiers.end());

    std::vector<TangleView> tangles;
    std::vector<bool> used(sdg.nodes.size());
    std::map<sgNodeID_t , bool> fused;

    for (auto nnid=0;nnid<sdg.nodes.size();++nnid) {
        if (used[nnid]) continue;
        if (sdg.nodes[nnid].status==NodeStatus::Deleted) continue;
        for(auto nid:{nnid,-nnid}){
            if (used[llabs(nid)] or fused[nid]) continue;
            TangleView t(this, {}, {});
            std::set<sgNodeID_t> to_explore = {nid};
            std::set<sgNodeID_t> explored;
            std::set<sgNodeID_t> frontier_ids;
            if (!fset.count(llabs(nid))) {
                //std::cout<<"New tangle from internal node "<<nid<<std::endl;
                t.internals.emplace_back(sdg.get_nodeview(nid));
                used[nid] = true;
            }
            else {
                //std::cout<<"New tangle from frontier node "<<nid<<std::endl;
                frontier_ids.emplace(nid);
                fused[nid]=true;
            }

            while (to_explore.size() > 0) {
                //std::cout<<" exploring "<<to_explore.size()<<" nodes"<<std::endl;
                std::set<sgNodeID_t> new_to_explore;
                for (auto &inid:to_explore) {
                    bool f = fset.count(llabs(inid));
                    //std::cout<<"  now exploring "<<inid<<(f?" as frontier":" as internal")<<std::endl;
                    for (auto &l:links[llabs(inid)]) {
                        if (f and l.source != inid) continue; //frontiers only explore on same end
                        if (explored.count(l.dest) > 0) continue;
                        if (fset.count(llabs(l.dest))) {
                            frontier_ids.emplace(l.dest);
                            fused[l.dest]=true;
                        } else {
                            if (used[llabs(l.dest)]) continue;
                            t.internals.emplace_back(this, llabs(l.dest));
                            used[llabs(l.dest)] = true;
                        }
                        if (explored.count(l.dest) == 0) new_to_explore.emplace(l.dest);
                    }
                    explored.emplace(inid);
                }
                std::swap(new_to_explore, to_explore);
            }
            for (auto &fid:frontier_ids) t.frontiers.emplace_back(sdg.get_nodeview(fid));
            if (include_disconnected or t.frontiers.size() > 0) tangles.emplace_back(t);
        }
        used[nnid] = true;
    }
    return tangles;
}

std::vector<TangleView> DistanceGraph::get_all_tangleviews(int fsize, float fmin_kci, float fmax_kci, bool include_disconnected) const {
//    std::cout<<"Tangles of size < "<<fsize<<" with "<<fmin_kci<<" < KCI < "<<fmax_kci<<" disconnected: "<<include_disconnected<<std::endl;
    std::vector<TangleView> tangles;
    std::vector<bool> used(sdg.nodes.size());
    for (auto nid=0;nid<sdg.nodes.size();++nid) {
        if (used[nid]) continue;
        auto &n=sdg.nodes[nid];
        if (n.status==NodeStatus::Deleted or n.sequence.size()>=fsize) continue;
        auto nv=get_nodeview(nid);
        if (fmin_kci>=0 and nv.kci()<fmin_kci) continue;
        if (fmax_kci>=0 and nv.kci()>fmax_kci) continue;
//        std::cout<<"New tangle from node "<<nid<<std::endl;
        TangleView t(this,{},{nid});
        used[nid]=true;
        std::set<sgNodeID_t> to_explore={nid};
        std::set<sgNodeID_t> explored;
        std::set<sgNodeID_t> frontier_ids;
        while(to_explore.size() > 0 ){
//            std::cout<<"  Exploring "<<to_explore.size()<<" nodes"<<std::endl;
            std::set<sgNodeID_t> new_to_explore;
            for (auto &inid:to_explore) {
//                std::cout<<"    Links from "<<inid<<std::endl;
                bool f=false;
                {
                    auto nv = get_nodeview(inid);//non-const for kci usage
                    if (nv.size() >= fsize and (fmin_kci < 0 or fmin_kci < nv.kci()) and
                        (fmax_kci < 0 or fmax_kci > nv.kci()))
                        f = true;
                }

                for (auto &l:links[llabs(inid)]){
                    if (f and l.source!=inid) continue; //frontiers only explore on same end
                    if (explored.count(l.dest)>0) continue;
                    auto nv = get_nodeview(l.dest);//non-const for kci usage
                    if (nv.size() >= fsize and (fmin_kci < 0 or fmin_kci < nv.kci()) and
                        (fmax_kci < 0 or fmax_kci > nv.kci())) {
                        frontier_ids.emplace(l.dest);
                    }
                    else {
                        if (used[llabs(l.dest)]) continue;
                        t.internals.emplace_back(this,llabs(l.dest));
                        used[llabs(l.dest)]=true;
                    }
                    if (explored.count(l.dest)==0) new_to_explore.emplace(l.dest);
                }

                explored.emplace(inid);
            }
            std::swap(new_to_explore,to_explore);
        }

        //Any frontier that is connected in both sides is actually internal
        std::set<sgNodeID_t> reclassified_frontiers;
        std::vector<sgNodeID_t> frontiers(frontier_ids.begin(),frontier_ids.end());
        for (auto i=0;i+1<frontiers.size();++i) {
            for (auto j=i+1;j<frontiers.size();++j) {
                if (frontiers[i]==-frontiers[j]) reclassified_frontiers.insert(llabs(frontiers[i]));
            }
        }
        std::set<sgNodeID_t> used_frontiers;
        for (auto &fnid:frontiers) {
            if (reclassified_frontiers.count(llabs(fnid))==0) t.frontiers.emplace_back(this,fnid);
            else if (fnid>0 and used_frontiers.count(fnid)==0) {
                t.internals.emplace_back(this,fnid);
                used_frontiers.emplace(fnid);
            }
        }
        if (include_disconnected or t.frontiers.size()>0) tangles.emplace_back(t);
    }
    return tangles;
}

std::ostream &operator<<(std::ostream &os, const DistanceGraph &dg) {
    uint64_t linkcount=0;
    for (auto lv:dg.links)
        for (auto l:lv)
            if (llabs(l.source) <= llabs(l.dest)) ++linkcount;

    os << "DistanceGraph "<< (dg.name.empty() ? "unnamed" : dg.name )<<" ("<<dg.sdg.name<<"): "<<linkcount<<" links";
    return os;
}

std::vector<uint64_t> DistanceGraph::nstats(std::vector<uint64_t> sizes, std::vector<int> stops){
    if (sizes.empty()) return std::vector<uint64_t>(stops.size(),0);
    uint64_t t=0;
    for (auto& s: sizes){
        t+=s;
    }

    uint64_t p=0;
    uint64_t next_stop=0;
    std::vector<uint64_t> nxx;
    // reverse sort sizes
    std::sort(sizes.begin(), sizes.end(), std::greater<>());
    for (auto& x: sizes){
        p+=x;
        while(next_stop<=stops.size() and p>=t*stops[next_stop]/100){
            nxx.push_back(x);
            next_stop+=1;
        }
    }
    return nxx;
}

std::string DistanceGraph::stats_by_kci() {

//    // Check the kci peak value is not -1
    if (sdg.ws.kmer_counters[0].get_kci_peak() < 0){
        return "KCI peak not set!";
    }

    std::vector<uint64_t> nokci_sizes;
    std::vector<uint64_t> kci0_sizes;
    std::vector<uint64_t> kci1_sizes;
    std::vector<uint64_t> kci2_sizes;
    std::vector<uint64_t> kci3_sizes;
    std::vector<uint64_t> kci4_sizes;

    uint64_t nokci_tips=0;
    uint64_t kci0_tips=0;
    uint64_t kci1_tips=0;
    uint64_t kci2_tips=0;
    uint64_t kci3_tips=0;
    uint64_t kci4_tips=0;

    uint64_t nokci_cr=0;
    uint64_t kci0_cr=0;
    uint64_t kci1_cr=0;
    uint64_t kci2_cr=0;
    uint64_t kci3_cr=0;
    uint64_t kci4_cr=0;

    uint64_t nokci_bubs=0;
    uint64_t kci0_bubs=0;
    uint64_t kci1_bubs=0;
    uint64_t kci2_bubs=0;
    uint64_t kci3_bubs=0;
    uint64_t kci4_bubs=0;

    std::vector<uint64_t> binned_bps(61, 0);
    std::vector<uint64_t> binned_counts(61, 0);

    for (NodeView& nv: get_all_nodeviews()){
        float kci = nv.kci();

        if (kci == -1){
            nokci_sizes.push_back(nv.size());
            if (nv.is_tip()) nokci_tips++;
            if (nv.is_canonical_repeat()) nokci_cr++;
            if (nv.is_bubble_side()) nokci_bubs++;
        } else if (kci<.5){
            kci0_sizes.push_back(nv.size());
            if (nv.is_tip()) kci0_tips++;
            if (nv.is_canonical_repeat()) kci0_cr++;
            if (nv.is_bubble_side()) kci0_bubs++;
        } else if (kci<=1.5){
            kci1_sizes.push_back(nv.size());
            if (nv.is_tip()) kci1_tips++;
            if (nv.is_canonical_repeat()) kci1_cr++;
            if (nv.is_bubble_side()) kci1_bubs++;
        } else if (kci<=2.5){
            kci2_sizes.push_back(nv.size());
            if (nv.is_tip()) kci2_tips++;
            if (nv.is_canonical_repeat()) kci2_cr++;
            if (nv.is_bubble_side()) kci2_bubs++;
        } else if (kci<=3.5){
            kci3_sizes.push_back(nv.size());
            if (nv.is_tip()) kci3_tips++;
            if (nv.is_canonical_repeat()) kci3_cr++;
            if (nv.is_bubble_side()) kci3_bubs++;
        } else {
            kci4_sizes.push_back(nv.size());
            if (nv.is_tip()) kci4_tips++;
            if (nv.is_canonical_repeat()) kci4_cr++;
            if (nv.is_bubble_side()) kci4_bubs++;
        }

        if (kci != -1) binned_bps[std::min((int)(kci*10),60)]+=nv.size();
        if (kci != -1) binned_counts[std::min((int)(kci*10),60)]++;

    }
    std::vector<uint64_t> all_sizes;
    all_sizes.insert( all_sizes.end(), nokci_sizes.begin(), nokci_sizes.end() );
    all_sizes.insert( all_sizes.end(), kci0_sizes.begin(), kci0_sizes.end() );
    all_sizes.insert( all_sizes.end(), kci1_sizes.begin(), kci1_sizes.end() );
    all_sizes.insert( all_sizes.end(), kci2_sizes.begin(), kci2_sizes.end() );
    all_sizes.insert( all_sizes.end(), kci3_sizes.begin(), kci3_sizes.end() );
    all_sizes.insert( all_sizes.end(), kci4_sizes.begin(), kci4_sizes.end() );

    uint64_t all_tips = nokci_tips+kci0_tips+kci1_tips+kci2_tips+kci3_tips+kci4_tips;
    uint64_t all_cr = nokci_cr+kci0_cr+kci1_cr+kci2_cr+kci3_cr+kci4_cr;
    uint64_t all_bubs = nokci_bubs+kci0_bubs+kci1_bubs+kci2_bubs+kci3_bubs+kci4_bubs;

    auto nstats_nokci = nstats(nokci_sizes);
    auto nstats_0 = nstats(kci0_sizes);
    auto nstats_1 = nstats(kci1_sizes);
    auto nstats_2 = nstats(kci2_sizes);
    auto nstats_3 = nstats(kci3_sizes);
    auto nstats_4 = nstats(kci4_sizes);
    auto all_stats = nstats(all_sizes);
    char buffer[10000]={0};
    std::sprintf(buffer," -------------------------------------------------------------------------------------------------------\n");
    std::sprintf(buffer+std::strlen(buffer),"|  KCI  |    Total bp   |   Nodes  |  Tips   | C Reps  | B Sides |     N25    |     N50    |     N75    |\n");
    std::sprintf(buffer+std::strlen(buffer),"|-------+---------------+----------+---------+---------+---------+------------+------------+------------|\n");
    std::sprintf(buffer+std::strlen(buffer),"| None  | %13lld | %8lld | %7lld | %7lld | %7lld | %10lld | %10lld | %10lld |\n", std::accumulate(nokci_sizes.begin(), nokci_sizes.end(), (uint64_t) 0), nokci_sizes.size(), nokci_tips, nokci_cr, nokci_bubs, nstats_nokci[0], nstats_nokci[1], nstats_nokci[2]);
    std::sprintf(buffer+std::strlen(buffer),"| < 0.5 | %13lld | %8lld | %7lld | %7lld | %7lld | %10lld | %10lld | %10lld |\n", std::accumulate(kci0_sizes.begin(), kci0_sizes.end(), (uint64_t) 0), kci0_sizes.size(), kci0_tips, kci0_cr, kci0_bubs, nstats_0[0], nstats_0[1], nstats_0[2]);
    std::sprintf(buffer+std::strlen(buffer),"| ~ 1   | %13lld | %8lld | %7lld | %7lld | %7lld | %10lld | %10lld | %10lld |\n", std::accumulate(kci1_sizes.begin(), kci1_sizes.end(), (uint64_t) 0), kci1_sizes.size(), kci1_tips, kci1_cr, kci1_bubs,  nstats_1[0], nstats_1[1], nstats_1[2]);
    std::sprintf(buffer+std::strlen(buffer),"| ~ 2   | %13lld | %8lld | %7lld | %7lld | %7lld | %10lld | %10lld | %10lld |\n", std::accumulate(kci2_sizes.begin(), kci2_sizes.end(), (uint64_t) 0), kci2_sizes.size(), kci2_tips, kci2_cr, kci2_bubs,  nstats_2[0], nstats_2[1], nstats_2[2]);
    std::sprintf(buffer+std::strlen(buffer),"| ~ 3   | %13lld | %8lld | %7lld | %7lld | %7lld | %10lld | %10lld | %10lld |\n", std::accumulate(kci3_sizes.begin(), kci3_sizes.end(), (uint64_t) 0), kci3_sizes.size(), kci3_tips, kci3_cr, kci3_bubs,  nstats_3[0], nstats_3[1], nstats_3[2]);
    std::sprintf(buffer+std::strlen(buffer),"| > 3.5 | %13lld | %8lld | %7lld | %7lld | %7lld | %10lld | %10lld | %10lld |\n", std::accumulate(kci4_sizes.begin(), kci4_sizes.end(), (uint64_t) 0), kci4_sizes.size(), kci4_tips, kci4_cr, kci4_bubs,  nstats_4[0], nstats_4[1], nstats_4[2]);
    std::sprintf(buffer+std::strlen(buffer),"|-------+---------------+----------+---------+---------+---------+------------+------------+------------|\n");
    std::sprintf(buffer+std::strlen(buffer),"| All   | %13lld | %8lld | %7lld | %7lld | %7lld | %10lld | %10lld | %10lld |\n", std::accumulate(all_sizes.begin(), all_sizes.end(), (uint64_t) 0), all_sizes.size(), all_tips, all_cr, all_bubs, all_stats[0], all_stats[1], all_stats[2]);
    std::sprintf(buffer+std::strlen(buffer)," -------------------------------------------------------------------------------------------------------\n");
    return std::string(buffer);
}

std::string DistanceGraph::simple_structure_stats() const{
    char buffer[10000]={0};
    uint64_t total=0;
    uint64_t canreps=0;
    uint64_t bubsides=0;
    uint64_t tips=0;

    for (const auto& nv: get_all_nodeviews()){
        total++;
        if (nv.is_canonical_repeat()) canreps++;
        if (nv.is_bubble_side()) bubsides++;
        if (nv.is_tip()) tips++;
    }
    std::sprintf(buffer, "%lld tips, %lld canonical repeats and %lld bubble sides out of %lld nodes\n", tips, canreps, bubsides, total);
    return std::string(buffer);
}

std::unordered_set<sgNodeID_t> DistanceGraph::get_connected_component(sgNodeID_t nid, bool signed_nodes) const{
    std::unordered_set<sgNodeID_t > component_nodes = {nid};
    std::unordered_set<sgNodeID_t > last_added = {nid};

    // while nodes are still beaing added
    while (last_added.size()>0){
        std::unordered_set<sgNodeID_t > new_last_added;
        // For all the outer layer of discovered nodes
        for (const auto& node: last_added){
            auto nv = get_nodeview(node);
            // Check FW and BW for new nodes and add to the collection and to the outer layer list
            for (const auto& o: nv.next()){
                auto onid = o.node().node_id();
                if (!signed_nodes) onid = abs(onid);
                if (component_nodes.find(onid)==component_nodes.end()){
                    component_nodes.insert(onid);
                    new_last_added.insert(onid);
                }
            }
            for (const auto& o: nv.prev()){
                auto onid = o.node().node_id();
                if (!signed_nodes) onid = abs(onid);
                if (component_nodes.find(onid)==component_nodes.end()){
                    component_nodes.insert(onid);
                    new_last_added.insert(onid);
                }
            }
        }
        last_added = new_last_added;
    }
    return component_nodes;
}

std::vector<std::unordered_set<sgNodeID_t>> DistanceGraph::get_all_connected_components(uint64_t min_nodes) const {
    std::vector<std::unordered_set<sgNodeID_t>> ccs;
    std::vector<bool> used(sdg.nodes.size());
    for (auto &nv:get_all_nodeviews()){
        if (used[nv.node_id()]) continue;
        ccs.push_back(get_connected_component(nv.node_id(),false));
        for (auto const &nid: ccs.back()) used[nid]=true;
        if (ccs.back().size()<min_nodes) ccs.pop_back();
    }
    return ccs;
}