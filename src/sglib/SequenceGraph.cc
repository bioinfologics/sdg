//
// Created by Bernardo Clavijo (EI) on 18/10/2017.
//

#include <iostream>
#include <fstream>
#include <sstream>
#include <set>
#include <math.h>
#include <sglib/logger/OutputLog.h>
#include "SequenceGraph.hpp"
#include "sglib/factories/KMerFactory.h"

bool Node::is_canonical() {
    for (size_t i=0,j=sequence.size()-1;i<j;++i,--j){
        char f=sequence[i];
        char r=sequence[j];
        switch(r){
            case 'A':
                r='T';
                break;
            case 'C':
                r='G';
                break;
            case 'G':
                r='C';
                break;
            case 'T':
                r='A';
                break;
            case 'N':
                break;
            default:
            std::cout<<"unexpected character in fasta file: '"<<r<<"'"<<std::endl;
        }
        if (f<r) return true;
        if (r<f) return false;
    }
    return true;
};

void Node::make_rc() {
    std::string rseq;
    rseq.resize(sequence.size());
    for (size_t i=0,j=sequence.size()-1;i<sequence.size();++i,--j){
        switch(sequence[j]){
            case 'A':
                rseq[i]='T';
                break;
            case 'C':
                rseq[i]='G';
                break;
            case 'G':
                rseq[i]='C';
                break;
            case 'T':
                rseq[i]='A';
                break;
        }
    }
    std::swap(sequence,rseq);
};

bool SequenceGraph::is_sane() {
    for (auto n=0;n<nodes.size();++n){
        for (auto l:links[n]){
            bool found=false;
            for (auto &dl:links[llabs(l.dest)]) {
                if (dl.dest == l.source and dl.source == l.dest){
                    found=true;
                    break;
                }
            }
            if (!found) return false;
        }
    }
}

sgNodeID_t SequenceGraph::add_node(Node n) {
    nodes.emplace_back(n);
    links.emplace_back();
    return (sgNodeID_t) nodes.size()-1;
}

void SequenceGraph::remove_node(sgNodeID_t n) {
    sgNodeID_t node=(n>0? n:-n);
    auto oldlinks=links[node];//this creates a copy to allow the iteration
    for (auto &l:oldlinks) remove_link(l.source,l.dest);
    nodes[node].status=sgNodeDeleted;
    //TODO: this is a lazy solution
    nodes[node].sequence.clear();
    //TODO: remove read mappings
}

void SequenceGraph::add_link(sgNodeID_t source, sgNodeID_t dest, int32_t d) {
    Link l(source,dest,d);
    links[(source > 0 ? source : -source)].emplace_back(l);
    std::swap(l.source,l.dest);
    links[(dest > 0 ? dest : -dest)].emplace_back(l);
}

Link SequenceGraph::get_link(sgNodeID_t source, sgNodeID_t dest) {
    for (auto l:links[(source > 0 ? source : -source)]) {
        if (l.source==source,l.dest==dest) return l;
    }
    return Link(0,0,0);
}

void SequenceGraph::remove_link(sgNodeID_t source, sgNodeID_t dest) {
    auto & slinks = links[(source > 0 ? source : -source)];
    slinks.erase(std::remove(slinks.begin(), slinks.end(), Link(source,dest,0)), slinks.end());
    auto & dlinks = links[(dest > 0 ? dest : -dest)];
    dlinks.erase(std::remove(dlinks.begin(), dlinks.end(), Link(dest,source,0)), dlinks.end());

}

std::vector<Link> SequenceGraph::get_fw_links( sgNodeID_t n){
    std::vector<Link> r;
    for (auto &l:links[(n>0 ? n : -n)]) if (l.source==-n) r.emplace_back(l);
    return r;
}

size_t SequenceGraph::count_active_nodes() {
    size_t t=0;
    for (auto &n: nodes){
        if (n.status==sgNodeStatus_t::sgNodeActive) ++t;
    }
    return t;
}

bool Link::operator==(const Link a){
    if (a.source == this->source && a.dest == this->dest){
        return true;
    }
    return false;
}

bool Link::operator<(const Link a) const {
    if (a.source < this->source){
        return true;
    } if (a.source == this->source && a.dest < this->dest) {
        return  true;
    }
    return false;
}

std::vector<std::vector<sgNodeID_t >> SequenceGraph::find_bubbles(std::vector<sgNodeID_t> component){
    std::vector<std::vector<sgNodeID_t >> bubbles;
    std::vector<sgNodeID_t > checked;
    // loop over all links in component, if 2 (or x) nodes have same source and dest, they are bubbles
    for (auto n: component){
        //std::cout << "n " << "name: " << oldnames[n] << std::endl;
        std::vector<sgNodeID_t > bubble;
        bubble.push_back(n);
        // if we haven't checked this node
        if (std::find(checked.begin(), checked.end(), n) == checked.end()){
            auto links_n = links[n];
            std::set<Link> links_set;
            std::vector<Link> links_uniq;

            for (auto link:links_n){
                auto s = link.source > 0 ? link.source:-link.source;
                auto d = link.dest > 0 ? link.dest:-link.dest;

                //std::cout << "source: " << oldnames[s] << " dest: " << oldnames[d] << std::endl;
                if (links_set.find(link) == links_set.end()){
                    //std::cout << "not duplicate " <<std::endl;
                    links_uniq.push_back(link);
                    links_set.insert(link);
                }
            }
            /*for (auto u:  links_uniq){
                //std::cout << oldnames[u.source] << " " << oldnames[u.dest] << " \n";
                std::cout << u.source << " " << u.dest << " \n";
            }
            std::cout << "\nlinks unique: " << links_uniq.size() << " links set " << links_set.size() << std::endl;*/
            // if l has exactly 2 links, one source and one dest, it may be a bubble
            if (links_uniq.size() == 2) {
                std::vector<sgNodeID_t> linked_to;
                std::map<sgNodeID_t, int> linked_2nd_degree;
                for (auto l:links_uniq){
                    if (l.dest == n or l.dest == -n){
                        auto s = l.source > 0 ? l.source:-l.source;
                        linked_to.push_back(s);
                        std::set<sgNodeID_t > seen_s2;
                        std::set<sgNodeID_t > seen_d;
                        for (auto l_2:links[s]){
                            auto s2 = l_2.source > 0 ? l_2.source:-l_2.source;
                            if (seen_s2.find(s2) == seen_s2.end()) {
                                if (s2 != s && s2 != n) {
                                    linked_2nd_degree[s2] += 1;
                                }
                                seen_s2.insert(s2);
                            }

                            auto d = l_2.dest > 0 ? l_2.dest:-l_2.dest;
                            if (seen_d.find(d) == seen_d.end()) {
                                if (d != s && d != n) {

                                    linked_2nd_degree[d] += 1;
                                }
                                seen_d.insert(d);
                            }
                        }
                    } else if (l.source == n or l.source == -n) {
                        auto s = l.dest > 0 ? l.dest:-l.dest;
                        linked_to.push_back(s);

                        std::set<sgNodeID_t > seen_s2;
                        std::set<sgNodeID_t > seen_d;
                        for (auto l_2:links[s]){
                            auto s2 = l_2.source > 0 ? l_2.source:-l_2.source;
                            if (seen_s2.find(s2) == seen_s2.end()) {

                                if (s2 !=s && s2 !=n){
                                linked_2nd_degree[s2] += 1;

                            }
                                seen_s2.insert(s2);
                            }

                             auto d = l_2.dest > 0 ? l_2.dest:-l_2.dest;
                            if (seen_d.find(d) == seen_d.end()) {
                                if (d != s && d != n) {
                                linked_2nd_degree[d] += 1;
                                //checked.push_back(d);
                            }
                                seen_d.insert(d);
                            }
                        }
                    }

                }

/*for (auto link2:linked_2nd_degree){

    auto s = link2.first > 0 ? link2.first:-link2.first;
    std::cout << "2nd degree l: " << oldnames[s] << " " << link2.second << " l first: " << link2.first << std::endl;
}*/
                    // if n is a bubble, other bubble contigs will share source and dest
                    // nope... but each bubble contig should occur twice...


                    // if an element is 2nd degree joined to n twice, and linked to the same nodes as n, its in a bubble
                    for (auto j:linked_2nd_degree){
                        if (j.second == 2){
                            checked.push_back(j.first);
                            auto links_j = links[j.first];
                            std::set<sgNodeID_t > joined_j;
                            for (auto l_j:links_j){
                                if (l_j.source == j.first or l_j.source == -j.first){
                                    auto l_j_abs = l_j.dest > 0? l_j.dest: -l_j.dest;
                                    joined_j.insert(l_j_abs);
                                } else if (l_j.dest == j.first or l_j.dest == -j.first){
                                    auto l_j_abs = l_j.source > 0? l_j.source: -l_j.source;

                                    joined_j.insert(l_j_abs);

                                }
                            }

                            auto s = j.first > 0 ? j.first:-j.first;
                            /*std::cout << "j: " << oldnames[s] << " id: " << j.first;
                            for (auto jhg:joined_j){
                                auto d = jhg > 0 ? jhg:-jhg;

                                std::cout << " "<< oldnames[d] << " ";
                            }
                            std::cout << std::endl;*/
                            if (joined_j.size() == linked_to.size()){
                                std::vector<sgNodeID_t > joined_vec;
                                for (auto joined:joined_j){
                                    joined_vec.push_back(joined);
                                    //std::cout << joined << " ";
                                }
                                //std::cout <<std::endl;
                                for (auto joined: linked_to){
                                    joined_vec.push_back(joined);
                                    //std::cout << joined << " ";
                                }
                                //std::cout << std::endl;

                                if (std::is_permutation(linked_to.begin(), linked_to.end(), joined_vec.begin()) ) {
                                    bubble.push_back(j.first);
                                }
                            }
                        }

                    }


            }

        }
        if (bubble.size() > 1){
            bubbles.push_back(bubble);
        }
    }
    std::cout << "Found " << bubbles.size() << " bubbles in component of " << component.size() << " nodes " <<std::endl;
    return bubbles;
}


std::vector<std::vector<sgNodeID_t>> SequenceGraph::connected_components(int max_nr_totalinks, int max_nr_dirlinks,
                                                                         int min_rsize) {
    std::vector<bool> used(nodes.size());
    std::vector<std::vector<sgNodeID_t>> components;
    //TODO: first find all repeats, add them as independent components and mark them as used.
    size_t max_component = 0;
    for (sgNodeID_t start_node=1;start_node<nodes.size();++start_node){
        if (false==used[start_node]){
            used[start_node]=true;
            //if start node is repeat, just add a single-node component
            std::set<sgNodeID_t> in_component;
            in_component.insert(start_node);
            std::set<sgNodeID_t> to_explore;
            to_explore.insert(start_node);
            while(to_explore.size()){
                auto n=*to_explore.begin();
                to_explore.erase(n);
                //find all connections that are not already used, include both in component and in to_explore
                for (auto l:links[n]) {
                    sgNodeID_t next=(l.dest>0 ? l.dest : -l.dest);
                    if (in_component.count(next)==0) {
                        in_component.insert(next);
                        to_explore.insert(next);
                    }
                }
            }
            components.emplace_back();
            for (sgNodeID_t n:in_component){
                components.back().push_back(n);
                used[n]=true;
            }
            if (in_component.size() > max_component){
                max_component = in_component.size();
            }
        }
    }
    std::cout << "Max component size: " << max_component << std::endl;
    return components;
}

void SequenceGraph::write(std::ofstream & output_file) {
    uint64_t count;
    count=nodes.size();
    output_file.write((char *) &count,sizeof(count));
    for (auto &n:nodes){
        output_file.write((char *) &n.status,sizeof(n.status));
        count=n.sequence.size();
        output_file.write((char *) &count,sizeof(count));
        output_file.write((char *) n.sequence.data(),count);
    }
    count=links.size();
    output_file.write((char *) &count,sizeof(count));
    for (auto nl:links) {
        uint64_t lcount=nl.size();
        output_file.write((char *) &lcount,sizeof(lcount));
        output_file.write((char *) nl.data(), sizeof(Link) * lcount);
    }
}

void SequenceGraph::read(std::ifstream & input_file) {
    uint64_t count;
    input_file.read((char *) &count,sizeof(count));
    nodes.clear();
    nodes.reserve(count);
    uint64_t active=0;
    for (auto i=0;i<count;++i){
        uint64_t seqsize;
        std::string seq;
        sgNodeStatus_t status;
        input_file.read((char *) &status,sizeof(status));
        input_file.read((char *) &seqsize,sizeof(seqsize));
        seq.resize(seqsize);
        input_file.read((char *) seq.data(),seqsize);
        nodes.emplace_back(seq);
        nodes.back().status=status;
        if (nodes.back().status==sgNodeStatus_t::sgNodeActive) ++active;
    }
    input_file.read((char *) &count,sizeof(count));
    links.resize(count);
    for (auto &nl:links) {
        uint64_t lcount;
        input_file.read((char *) &lcount,sizeof(lcount));
        nl.resize(lcount,{0,0,0});
        input_file.read((char *) nl.data(), sizeof(Link) * lcount);
    }
}

void SequenceGraph::load_from_gfa(std::string filename) {
    std::string line;
    this->filename=filename;
    //check the filename ends in .gfa
    if (filename.size()>4 and filename.substr(filename.size()-4,4)==".gfa"){
        this->fasta_filename=filename.substr(0,filename.size()-4)+".fasta";
    }
    else throw std::invalid_argument("filename of the gfa input does not end in gfa, it ends in '"+filename.substr(filename.size()-4,4)+"'");

    std::ifstream gfaf(filename);
    if (!gfaf) throw std::invalid_argument("Can't read gfa file");
    std::getline(gfaf, line);
    if (line!="H\tVN:Z:1.0") std::cout<<"WARNING, first line of gfa doesn't correspond to GFA1"<<std::endl;

    std::ifstream fastaf(fasta_filename);
    std::cout << "fasta filesname: " << fasta_filename << std::endl;
    if (!fastaf) throw std::invalid_argument("Can't read fasta file");



    //load all sequences from fasta file if they're not canonical, flip and remember they're flipped
    std::cout<<"Loading sequences from "<<fasta_filename<<std::endl;

    std::string name,seq="";
    seq.reserve(10000000); //stupid hack but probably useful to reserve
    oldnames_to_ids.clear();
    oldnames.push_back("");
    nodes.clear();
    links.clear();
    add_node(Node("",sgNodeStatus_t::sgNodeDeleted)); //an empty deleted node on 0, just to skip the space
    sgNodeID_t nextid=1;
    uint64_t rcnodes=0;
    while(!fastaf.eof()){
        std::getline(fastaf,line);
        if (fastaf.eof() or line[0]=='>'){
            if (!name.empty()) {
                //rough ansi C and C++ mix but it works
                if (oldnames_to_ids.find(name) != oldnames_to_ids.end())
                    throw std::logic_error("sequence " + name + " is already defined");
                oldnames_to_ids[name] = add_node(Node(seq));
                oldnames.push_back(name);
                //reverse seq if not canonical
                if (!nodes.back().is_canonical()) {
                    nodes.back().make_rc();
                    oldnames_to_ids[name] = -oldnames_to_ids[name];
                    ++rcnodes;
                }
            }
            name.clear();
            for (auto i = 1; i < line.size() and line[i] != ' '; ++i) name += line[i];
            seq = "";
        } else {
            seq+=line;
        }
    }
    std::cout<<nodes.size()-1<<" nodes loaded! "<<rcnodes<<" canonised"<<std::endl;

    //load store all conections.

    std::string gfa_rtype,gfa_source,gfa_sourcedir,gfa_dest,gfa_destdir,gfa_cigar;
    sgNodeID_t src_id,dest_id;
    int32_t dist;
    uint64_t lcount=0;
    uint64_t dist_egt0(0);
    while(std::getline(gfaf, line) and !gfaf.eof()) {
        std::istringstream iss(line);
        iss >> gfa_rtype;
        if (gfa_rtype == "L"){
            iss >> gfa_source;
            iss >> gfa_sourcedir;
            iss >> gfa_dest;
            iss >> gfa_destdir;
            iss >> gfa_cigar;
            //std::cout<<"'"<<source<<"' '"<<gfa_sourcedir<<"' '"<<dest<<"' '"<<destdir<<"'"<<std::endl;
            if (oldnames_to_ids.find(gfa_source)==oldnames_to_ids.end()){
                oldnames_to_ids[gfa_source] = add_node(Node(""));
                //std::cout<<"added source!" <<source<<std::endl;
            }
            if (oldnames_to_ids.find(gfa_dest)==oldnames_to_ids.end()){
                oldnames_to_ids[gfa_dest] = add_node(Node(""));
                //std::cout<<"added dest! "<<dest<<std::endl;
            }
            src_id=oldnames_to_ids[gfa_source];
            dest_id=oldnames_to_ids[gfa_dest];

            //Go from GFA's "Links as paths" to a normal "nodes have sinks (-/start/left) and sources (+/end/right)
            if (gfa_sourcedir=="+") src_id=-src_id;
            if (gfa_destdir=="-") dest_id=-dest_id;
            dist=0;
            if (gfa_cigar.size()>1 and gfa_cigar[gfa_cigar.size()-1]=='M') {
                //TODO: better checks for M
                dist=-atoi(gfa_cigar.c_str());
                if (dist>=0) {
                    dist_egt0++;
                }
            }
            add_link(src_id,dest_id,dist);
            ++lcount;
        }
    }
    if (dist_egt0 > lcount*0.5f) {
        sglib::OutputLog() << "Warning: The loaded graph contains " << dist_egt0 << " non-overlapping links out of " << lcount << std::endl;
    }
    sglib::OutputLog() <<nodes.size()-1<<" nodes after connecting with "<<lcount<<" links"<<std::endl;
}

void SequenceGraph::write_to_gfa(std::string filename, const std::unordered_set<sgNodeID_t> & mark_red, const std::vector<double> & depths,
                                 const std::unordered_set<sgNodeID_t> & selected_nodes={}, const std::vector<std::vector<Link>> & arg_links){
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
    std::cout<<"Writing sequences to "<<fasta_filename<<std::endl;

    for (sgNodeID_t i=1;i<nodes.size();++i){
        if (nodes[i].status==sgNodeDeleted) continue;
        if (!selected_nodes.empty() and selected_nodes.count(i)==0 and selected_nodes.count(-i)==0) continue;
        fastaf<<">seq"<<i<<std::endl<<nodes[i].sequence<<std::endl;
        gfaf<<"S\tseq"<<i<<"\t*\tLN:i:"<<nodes[i].sequence.size()<<"\tUR:Z:"<<fasta_filename
                <<(mark_red.count(i)?"\tCL:Z:red":"")<<(depths.empty() or isnan(depths[i])?"":"\tDP:f:"+std::to_string(depths[i]))<<std::endl;
    }

    for (auto &ls:(arg_links.size()>0? arg_links:links)){
        for (auto &l:ls)
            if (l.source<=l.dest and (selected_nodes.empty() or
                    selected_nodes.count(l.source)>0 or selected_nodes.count(-l.source)>0 or
                            selected_nodes.count(l.dest)>0 or selected_nodes.count(-l.dest)>0)) {
                gfaf<<"L\t";
                if (l.source>0) gfaf<<"seq"<<l.source<<"\t-\t";
                else gfaf<<"seq"<<-l.source<<"\t+\t";
                if (l.dest>0) gfaf<<"seq"<<l.dest<<"\t+\t";
                else gfaf<<"seq"<<-l.dest<<"\t-\t";
                gfaf<<(l.dist<0 ? -l.dist : 0)<<"M"<<std::endl;
            }
    }

}

std::vector<sgNodeID_t> SequenceGraph::oldnames_to_nodes(std::string _oldnames) {
    std::vector<sgNodeID_t> nv;
    const char * s = _oldnames.c_str();
    const char * sign;
    std::string oldname;
    sgNodeID_t node;
    while (*s!= NULL){
        for (sign=s;*sign!=NULL and *sign !='+' and *sign!='-';++sign);
        if (sign ==s or *sign==NULL ) throw std::invalid_argument("invalid path specification");
        oldname=s;
        oldname.resize(sign-s);
        node=oldnames_to_ids[oldname];
        if (0==node) throw std::invalid_argument("node "+oldname+" doesn't exist in graph");
        if (*sign=='-') node=-node;
        nv.push_back(node);
        s=sign+1;
        if (*s==',') ++s;
        else if (*s!=NULL) throw std::invalid_argument("invalid path specification");

    }
    return nv;
}

std::string SequenceGraphPath::get_fasta_header() {
    std::string h=">sgPath_";
    for (auto &n:nodes) {
        h += std::to_string(n)+",";
    }
    h.resize(h.size()-1);
    return h;
}

std::string SequenceGraphPath::get_sequence() {
    std::string s="";
    sgNodeID_t pnode=0;
    // just iterate over every node in path - contig names are converted to ids at construction
    for (auto &n:nodes) {
        std::string nseq;
        if (n>0){
            nseq=sg.nodes[n].sequence;
        } else {
            auto rcn=sg.nodes[-n];
            rcn.make_rc();
            nseq=rcn.sequence;
        }
        if (pnode !=0){
            //find link between pnode' output (+pnode) and n's sink (-n)
            auto l=sg.links[(pnode>0 ? pnode:-pnode)].begin();
            for (;l!=sg.links[(pnode>0 ? pnode:-pnode)].end();++l)
                if (l->source==pnode and l->dest==n) break;
            if (l==sg.links[(pnode>0 ? pnode:-pnode)].end()) {
                std::cout<<"can't find a link between "<<pnode<<" and "<<n<<std::endl;
                throw std::runtime_error("path has no link");
            } else {
                if (l->dist>0){
                    for (auto c=l->dist;c>0;--c) s+="N";
                }
                else {
                    auto ovl=-l->dist;
                    for (auto s1=s.c_str()+s.size()-ovl,s2=nseq.c_str();*s1!=NULL;++s1,++s2)
                        if (*s1!=*s2)
                            throw std::runtime_error("path overlap is invalid!");
                    nseq.erase(0,ovl);
                }
            }
        }
        s+=nseq;
        pnode=-n;
    }
    return s;
}

size_t SequenceGraphPath::get_sequence_size_fast() {
    size_t size=0;
    //std::string s="";
    sgNodeID_t pnode=0;
    // just iterate over every node in path - contig names are converted to ids at construction
    for (auto &n:nodes) {
        std::string nseq;
        size=sg.nodes[llabs(n)].sequence.size();
        if (pnode !=0){
            //find link between pnode' output (+pnode) and n's sink (-n)
            auto l=sg.links[llabs(pnode)].begin();
            for (;l!=sg.links[llabs(pnode)].end();++l)
                if (l->source==pnode and l->dest==n) break;
            if (l==sg.links[llabs(pnode)].end()) {
                std::cout<<"can't find a link between "<<pnode<<" and "<<n<<std::endl;
                throw std::runtime_error("path has no link");
            } else {
                size+=l->dist;
            }
        }
        pnode=-n;
    }
    return size;
}

std::vector<SequenceGraphPath> SequenceGraph::get_all_unitigs(uint16_t min_nodes) {
    std::vector<SequenceGraphPath> unitigs;
    std::vector<bool> used(nodes.size(),false);

    for (auto n=1;n<nodes.size();++n){
        if (used[n] or nodes[n].status==sgNodeDeleted) continue;
        used[n]=true;
        SequenceGraphPath path(*this,{n});

        //two passes: 0->fw, 1->bw, path is inverted twice, so still n is +
        for (auto pass=0; pass<2; ++pass) {
            //walk til a "non-unitig" junction
            for (auto fn = get_fw_links(path.nodes.back()); fn.size() == 1; fn = get_fw_links(path.nodes.back())) {
                if (!used[llabs(fn[0].dest)] and get_bw_links(fn[0].dest).size() == 1) {
                    path.nodes.emplace_back(fn[0].dest);
                    used[fn[0].dest > 0 ? fn[0].dest : -fn[0].dest] = true;
                } else break;
            }
            path.reverse();
        }
        if (path.nodes.size()>=min_nodes) unitigs.push_back(path);
    }
    return unitigs;
}

uint32_t SequenceGraph::join_all_unitigs() {
    uint32_t joined=0;
    for (auto p:get_all_unitigs(2)){
        join_path(p);
        ++joined;
    }
    return joined;
}

void SequenceGraph::join_path(SequenceGraphPath p, bool consume_nodes) {
    std::set<sgNodeID_t> pnodes;
    for (auto n:p.nodes) {
        pnodes.insert( n );
        pnodes.insert( -n );
    }
    if (!p.is_canonical()) p.reverse();
    sgNodeID_t new_node=add_node(Node(p.get_sequence()));
    //TODO:check, this may have a problem with a circle
    for (auto l:get_bw_links(p.nodes.front())) add_link(new_node,l.dest,l.dist);

    for (auto l:get_fw_links(p.nodes.back())) add_link(-new_node,l.dest,l.dist);

    //TODO: update read mappings
    if (consume_nodes) {
        for (auto n:p.nodes) {
            //check if the node has neighbours not included in the path.
            bool ext_neigh=false;
            if (n!=p.nodes.back()) for (auto l:get_fw_links(n)) if (pnodes.count(l.dest)==0) ext_neigh=true;
            if (n!=p.nodes.front()) for (auto l:get_bw_links(n)) if (pnodes.count(l.dest)==0) ext_neigh=true;
            if (ext_neigh) continue;
            remove_node(n);
        }
    }
}

void SequenceGraphPath::reverse(){
    std::vector<sgNodeID_t> newn;
    for (auto n=nodes.rbegin();n<nodes.rend();++n) newn.emplace_back(-*n);
    //std::swap(nodes,newn);
    nodes=newn;
}

bool SequenceGraphPath::is_canonical() {
    auto rp=*this;
    rp.reverse();
    return this->get_sequence()<rp.get_sequence();
}

bool SequenceGraphPath::is_unitig() {
    for (auto i=0;i<nodes.size();++i) {
        auto fwl=sg.get_fw_links(nodes[i]);
        auto bwl=sg.get_bw_links(nodes[i]);
        if (i>0){
            if (bwl.size()!=1) return false;
            if (bwl[0].dest!=-nodes[i-1]) return false;
        }
        if (i<nodes.size()-1){
            if (fwl.size()!=1) return false;
            if (fwl[0].dest!=nodes[i+1]) return false;
        }
    }
    return true;
}

const bool SequenceGraphPath::operator<(const SequenceGraphPath &other) {
    for (auto i=0;i<nodes.size();++i){
        if (other.nodes.size()<i) return true;
        if (nodes[i]<other.nodes[i]) return true;
        if (nodes[i]>other.nodes[i]) return false;
    }
    return false;
}

const bool SequenceGraphPath::operator==(const SequenceGraphPath &other) {
    if (other.nodes.size()!=nodes.size()) return false;
    for (auto i=0;i<nodes.size();++i){
        if (nodes[i]!=other.nodes[i]) return false;
    }
    return true;
}

bool SequenceGraphPath::extend_if_coherent(SequenceGraphPath s) {
//    int offset=-1;
//    for (auto i=0;i<nodes.size();++i) {
//        if (nodes[i] == s.nodes[0]) {
//            offset = i;
//            break;
//        }
//    }
//    if (offset<=0) {
//    if (offset<0) {
//        offset=1;
//        for (auto i=0;i<nodes.size();++i) {
//            if (nodes[i] == s.nodes[0]) {
//                offset = -i;
//                break;
//            }
//        }
//    }
}

std::vector<SequenceSubGraph> SequenceGraph::get_all_tribbles() {


    for (sgNodeID_t n=1;n<nodes.size();++n) {
        //Heuristic to find "tribbles"
        // A --- B -- C -- H
        //  \     \      /
        //   \     E    /
        //    \     \  /
        //      F -- G
        // A->[B-F]
        // B->[C-E]
        // C->H
        // E->G
        // F->G
        // F->H
        auto a_fw=get_fw_links(n);
        if (a_fw.size()!=2) continue;
        auto b_fw=get_fw_links(a_fw[0].dist);

        sgNodeID_t A, B, C, D, E, F, G, H;
    }


}


void SequenceGraph::expand_node(sgNodeID_t nodeID, std::vector<std::vector<sgNodeID_t>> bw,
                                std::vector<std::vector<sgNodeID_t>> fw) {
    if (nodeID<0) {
        //std::cout<<"WARNING: ERROR: expand_node only accepts positive nodes!"<<std::endl;
        //return;
        nodeID=-nodeID;
        std::swap(bw,fw);
        for (auto &x:bw) for (auto &xx:x) xx=-xx;
        for (auto &x:fw) for (auto &xx:x) xx=-xx;
    }
    //TODO: check all inputs are included in bw and all outputs are included in fw, only once
    auto orig_links=links[nodeID];
    std::vector<std::vector<Link>> new_links;
    new_links.resize(bw.size());
    for (auto l:orig_links){
        for (auto in=0;in<bw.size();++in){
            if (l.source==nodeID and std::find(bw[in].begin(),bw[in].end(),l.dest)!=bw[in].end()) {
                new_links[in].push_back(l);
                break;
            }
        }
        for (auto in=0;in<fw.size();++in){
            if (l.source==-nodeID and std::find(fw[in].begin(),fw[in].end(),l.dest)!=fw[in].end()) {
                new_links[in].push_back(l);
                break;
            }
        }
    }
    //Create all extra copies of the node.
    for (auto in=0;in<new_links.size();++in){
        auto new_node=add_node(Node(nodes[llabs(nodeID)].sequence));
        for (auto l:new_links[in]) {
            add_link((l.source>0 ? new_node:-new_node),l.dest,l.dist);
        }
    }
    remove_node(nodeID);

}

std::vector<std::pair<sgNodeID_t,int64_t>> SequenceGraph::get_distances_to(sgNodeID_t n,
                                                                            std::set<sgNodeID_t> destinations,
                                                                            int64_t max_dist) {
    std::vector<std::pair<sgNodeID_t,int64_t>> current_nodes,next_nodes,final_nodes;
    current_nodes.push_back({n,-nodes[llabs(n)].sequence.size()});
    int rounds=20;
    while (not current_nodes.empty() and --rounds>0){
        //std::cout<<"starting round of context expansion, current nodes:  ";
        //for (auto cn:current_nodes) std::cout<<cn.first<<"("<<cn.second<<")  ";
        //std::cout<<std::endl;
        next_nodes.clear();
        for (auto nd:current_nodes){
            //if node is a destination add it to final nodes
            if (destinations.count(nd.first)>0 or destinations.count(-nd.first)>0){
                final_nodes.push_back(nd);
            }
            //else get fw links, compute distances for each, and if not >max_dist add to nodes
            else {
                for (auto l:get_fw_links(nd.first)){
                    int64_t dist=nd.second;
                    dist+=nodes[llabs(nd.first)].sequence.size();
                    dist+=l.dist;
                    //std::cout<<"candidate next node "<<l.dest<<" at distance "<<dist<<std::endl;
                    if (dist<=max_dist) next_nodes.push_back({l.dest,dist});
                }
            }
        }
        current_nodes=next_nodes;
    }
    return final_nodes;
}

std::vector<SequenceSubGraph> SequenceGraph::get_all_bubbly_subgraphs(uint32_t maxsubgraphs) {
    std::vector<SequenceSubGraph> subgraphs;
    std::vector<bool> used(nodes.size(),false);
    /*
     * the loop always keep the first and the last elements as c=2 collapsed nodes.
     * it starts with a c=2 node, and goes thorugh all bubbles fw, then reverts the subgraph and repeats
     */
    SequenceSubGraph subgraph(*this);
    for (auto n=1;n<nodes.size();++n){
        if (used[n] or nodes[n].status==sgNodeDeleted) continue;
        subgraph.nodes.clear();

        subgraph.nodes.push_back(n);
        bool circular=false;
        //two passes: 0->fw, 1->bw, path is inverted twice, so still n is +
        for (auto pass=0; pass<2; ++pass) {
            //while there's a possible bubble fw.
            for (auto fn = get_fw_links(subgraph.nodes.back()); fn.size() == 2; fn = get_fw_links(subgraph.nodes.back())) {
                //if it is not a real bubble, get out.
                if (get_bw_links(fn[0].dest).size()!=1 or get_bw_links(fn[0].dest).size()!=1) break;
                auto fl1=get_fw_links(fn[0].dest);
                if (fl1.size()!=1) break;
                auto fl2=get_fw_links(fn[1].dest);
                if (fl2.size()!=1) break;
                if (fl2[0].dest!=fl1[0].dest) break;
                auto next_end=fl2[0].dest;
                //all conditions met, update subgraph
                if (used[llabs(fn[0].dest)] or used[llabs(fn[1].dest)] or used[llabs(next_end)]) {
                    circular=true;
                    break;
                }
                subgraph.nodes.push_back(fn[0].dest);
                subgraph.nodes.push_back(fn[1].dest);
                subgraph.nodes.push_back(next_end);
                used[llabs(fn[0].dest)]=true;
                used[llabs(fn[1].dest)]=true;
                used[llabs(next_end)]=true;
                used[n]=true;
            }
            SequenceSubGraph new_subgraph(*this);
            for (auto it=subgraph.nodes.rbegin();it<subgraph.nodes.rend();++it) new_subgraph.nodes.push_back(-*it);
            std::swap(new_subgraph.nodes,subgraph.nodes);
        }
        if (subgraph.nodes.size()>6 and not circular) {
            subgraphs.push_back(subgraph);
            if (subgraphs.size()==maxsubgraphs) break;
            //std::cout<<"Bubbly path found: ";
            //for (auto &n:subgraph) std::cout<<"  "<<n<<" ("<<sg.nodes[(n>0?n:-n)].sequence.size()<<"bp)";
            //std::cout<<std::endl;
        }
    }
    return subgraphs;
}

void SequenceGraph::print_bubbly_subgraph_stats(const std::vector<SequenceSubGraph> &bubbly_paths) {
    std::vector<uint64_t> solved_sizes,original_sizes;
    sglib::OutputLog()<<bubbly_paths.size()<<" bubbly paths"<<std::endl;
    auto &log_no_date=sglib::OutputLog(sglib::LogLevels::INFO,false);
    uint64_t total_size=0,total_solved_size=0;
    for (auto &bp:bubbly_paths){
        total_size+=bp.total_size();
        SequenceGraphPath p1(*this),p2(*this);
        original_sizes.push_back(nodes[llabs(bp.nodes.front())].sequence.size());
        original_sizes.push_back(nodes[llabs(bp.nodes.back())].sequence.size());
        solved_sizes.push_back(nodes[llabs(bp.nodes.front())].sequence.size());
        total_solved_size+=solved_sizes.back();
        solved_sizes.push_back(nodes[llabs(bp.nodes.back())].sequence.size());
        total_solved_size+=solved_sizes.back();
        for (auto i=1; i<bp.nodes.size()-1; ++i){
            original_sizes.push_back(nodes[llabs(bp.nodes[i])].sequence.size());
            if (i%3==0){
                p1.nodes.push_back(bp.nodes[i]);
                p2.nodes.push_back(bp.nodes[i]);
            } else if (i%3==1) p1.nodes.push_back(bp.nodes[i]);
            else p2.nodes.push_back(bp.nodes[i]);
        }
        solved_sizes.push_back(p1.get_sequence().size());
        total_solved_size+=solved_sizes.back();
        solved_sizes.push_back(p2.get_sequence().size());
        total_solved_size+=solved_sizes.back();
    }
    std::sort(solved_sizes.rbegin(),solved_sizes.rend());
    std::sort(original_sizes.rbegin(),original_sizes.rend());
    auto on20s=total_size*.2;
    auto on50s=total_size*.5;
    auto on80s=total_size*.8;
    sglib::OutputLog() <<"Currently "<<original_sizes.size()<<" sequences with "<<total_size<<"bp, ";
    uint64_t acc=0;
    for (auto s:original_sizes){
        if (acc<on20s and acc+s>on20s) log_no_date<<"N20: "<<s<<"  ";
        if (acc<on50s and acc+s>on50s) log_no_date<<"N50: "<<s<<"  ";
        if (acc<on80s and acc+s>on80s) log_no_date<<"N80: "<<s<<"  ";
        acc+=s;
    }
    log_no_date<<std::endl;
    auto sn20s=total_solved_size*.2;
    auto sn50s=total_solved_size*.5;
    auto sn80s=total_solved_size*.8;
    sglib::OutputLog()<<"Potentially "<<solved_sizes.size()<<" sequences with "<<total_solved_size<<"bp, ";
    acc=0;
    for (auto s:solved_sizes){
        if (acc<sn20s and acc+s>sn20s) log_no_date<<"N20: "<<s<<"  ";
        if (acc<sn50s and acc+s>sn50s) log_no_date<<"N50: "<<s<<"  ";
        if (acc<sn80s and acc+s>sn80s) log_no_date<<"N80: "<<s<<"  ";
        acc+=s;
    }
    log_no_date<<std::endl;

}

const uint64_t SequenceSubGraph::total_size() const {
    uint64_t t=0;
    for (auto &n:nodes) t+=sg.nodes[llabs(n)].sequence.size();
    return t;
}

std::vector<SequenceGraphPath> SequenceGraph::find_all_paths_between(sgNodeID_t from,sgNodeID_t to, int64_t max_size, int max_nodes, bool abort_on_loops) {
    std::vector<SequenceGraphPath> current_paths,next_paths,final_paths;
    for(auto &fl:get_fw_links(from)) current_paths.emplace_back(SequenceGraphPath(*this,{fl.dest}));
    //int rounds=20;
    while (not current_paths.empty() and --max_nodes>0){
        //std::cout<<"starting round of context expansion, current nodes:  ";
        //for (auto cn:current_nodes) std::cout<<cn.first<<"("<<cn.second<<")  ";
        //std::cout<<std::endl;
        next_paths.clear();
        for (auto &p:current_paths){
            //if node is a destination add it to final nodes
            if (p.nodes.back()==-to) return {};//std::cout<<"WARNING: found path to -TO node"<<std::endl;
            if (p.nodes.back()==to) {
                final_paths.push_back(p);
                final_paths.back().nodes.pop_back();
            }
                //else get fw links, compute distances for each, and if not >max_dist add to nodes
            else {
                for (auto l:get_fw_links(p.nodes.back())){
                    if (abort_on_loops and std::find(p.nodes.begin(),p.nodes.end(),l.dest)!=p.nodes.end()) {
                        //std::cout<<"Loop detected, aborting pathing attempt!"<<std::endl;
                        return {};
                        //continue;
                    }
                    next_paths.push_back(p);
                    next_paths.back().nodes.push_back(l.dest);
                    if (l.dest!=to and next_paths.back().get_sequence_size_fast()>max_size) next_paths.pop_back();
                }
            }
        }
        current_paths=next_paths;
    }
    return final_paths;
}

void SequenceGraph::create_index() {
    kmer_to_graphposition.clear();
    class kmerPosFactory : protected KMerFactory {
    public:
        explicit kmerPosFactory(uint8_t k) : KMerFactory(k) {}

        ~kmerPosFactory() {
        }
        void setFileRecord(FastaRecord &rec) {
            currentRecord = rec;
            fkmer=0;
            rkmer=0;
            last_unknown=0;
        }

        // TODO: Adjust for when K is larger than what fits in uint64_t!
        const bool next_element(std::vector<std::pair<uint64_t,graphPosition>> &mers) {
            uint64_t p(0);
            graphPosition pos;
            while (p < currentRecord.seq.size()) {
                //fkmer: grows from the right (LSB)
                //rkmer: grows from the left (MSB)
                bases++;
                fillKBuf(currentRecord.seq[p], p, fkmer, rkmer, last_unknown);
                p++;
                if (last_unknown >= K) {
                    if (fkmer <= rkmer) {
                        // Is fwd
                        pos.node=currentRecord.id;
                        pos.pos=p;
                        mers.emplace_back(fkmer, pos);
                    } else {
                        // Is bwd
                        pos.node=-currentRecord.id;
                        pos.pos=p;
                        mers.emplace_back(rkmer, pos);
                    }
                }
            }
            return false;
        }

    private:
        FastaRecord currentRecord;
        uint64_t bases;
    };

    sglib::OutputLog(sglib::INFO) << "Indexing graph..."<<std::endl;
    const int k = 31;
    uint64_t total_k=0;
    std::vector<std::pair<uint64_t,graphPosition>> kidxv;
    for (auto &n:nodes) if (n.sequence.size()>=k) total_k+=n.sequence.size()+1-k;
    kidxv.reserve(total_k);
    FastaRecord r;
    kmerPosFactory kcf({k});
    for (sgNodeID_t n=1;n<nodes.size();++n){
        if (nodes[n].sequence.size()>=k){
            r.id=n;
            r.seq=nodes[n].sequence;
            kcf.setFileRecord(r);
            kcf.next_element(kidxv);
        }
    }
    sglib::OutputLog(sglib::INFO)<<kidxv.size()<<" kmers in total"<<std::endl;
    sglib::OutputLog(sglib::INFO) << "  Sorting..."<<std::endl;
    std::sort(kidxv.begin(),kidxv.end(),[](const std::pair<uint64_t,graphPosition> & a, const std::pair<uint64_t,graphPosition> & b){return a.first<b.first;});
    sglib::OutputLog(sglib::INFO) << "  Merging..."<<std::endl;
    auto wi=kidxv.begin();
    auto ri=kidxv.begin();
    auto nri=kidxv.begin();
    while (ri<kidxv.end()){
        while (nri!=kidxv.end() and nri->first==ri->first) ++nri;
        if (nri-ri==1) {
            *wi=*ri;
            ++wi;
        }
        ri=nri;
    }

    kidxv.resize(wi-kidxv.begin());
    sglib::OutputLog(sglib::INFO)<<kidxv.size()<<" unique kmers in index, creating map"<<std::endl;
    kmer_to_graphposition.insert(kidxv.begin(),kidxv.end());
    sglib::OutputLog(sglib::INFO)<<"map ready"<<std::endl;
}

void SequenceGraph::create_63mer_index() {
    k63mer_to_graphposition.clear();
    class kmerPosFactory128 : protected KMerFactory128 {
    public:
        explicit kmerPosFactory128(uint8_t k) : KMerFactory128(k) {}

        ~kmerPosFactory128() {
        }
        void setFileRecord(FastaRecord &rec) {
            currentRecord = rec;
            fkmer=0;
            rkmer=0;
            last_unknown=0;
        }

        // TODO: Adjust for when K is larger than what fits in uint64_t!
        const bool next_element(std::vector<std::pair<__uint128_t,graphPosition>> &mers) {
            uint64_t p(0);
            graphPosition pos;
            while (p < currentRecord.seq.size()) {
                //fkmer: grows from the right (LSB)
                //rkmer: grows from the left (MSB)
                bases++;
                fillKBuf(currentRecord.seq[p], p, fkmer, rkmer, last_unknown);
                p++;
                if (last_unknown >= K) {
                    if (fkmer <= rkmer) {
                        // Is fwd
                        pos.node=currentRecord.id;
                        pos.pos=p;
                        mers.emplace_back(fkmer, pos);
                    } else {
                        // Is bwd
                        pos.node=-currentRecord.id;
                        pos.pos=p;
                        mers.emplace_back(rkmer, pos);
                    }
                }
            }
            return false;
        }

    private:
        FastaRecord currentRecord;
        uint64_t bases;
    };

    sglib::OutputLog(sglib::INFO) << "Indexing graph..."<<std::endl;
    const int k = 63;
    uint64_t total_k=0;
    std::vector<std::pair<__uint128_t,graphPosition>> kidxv;
    for (auto &n:nodes) if (n.sequence.size()>=k) total_k+=n.sequence.size()+1-k;
    kidxv.reserve(total_k);
    FastaRecord r;
    kmerPosFactory128 kcf({k});
    for (sgNodeID_t n=1;n<nodes.size();++n){
        if (nodes[n].sequence.size()>=k){
            r.id=n;
            r.seq=nodes[n].sequence;
            kcf.setFileRecord(r);
            kcf.next_element(kidxv);
        }
    }
    sglib::OutputLog(sglib::INFO)<<kidxv.size()<<" kmers in total"<<std::endl;
    sglib::OutputLog(sglib::INFO) << "  Sorting..."<<std::endl;
    std::sort(kidxv.begin(),kidxv.end(),[](const std::pair<__uint128_t,graphPosition> & a, const std::pair<__uint128_t,graphPosition> & b){return a.first<b.first;});
    sglib::OutputLog(sglib::INFO) << "  Merging..."<<std::endl;
    auto wi=kidxv.begin();
    auto ri=kidxv.begin();
    auto nri=kidxv.begin();
    while (ri<kidxv.end()){
        while (nri!=kidxv.end() and nri->first==ri->first) ++nri;
        if (nri-ri==1) {
            *wi=*ri;
            ++wi;
        }
        ri=nri;
    }

    kidxv.resize(wi-kidxv.begin());
    sglib::OutputLog(sglib::INFO)<<kidxv.size()<<" unique kmers in index, creating map"<<std::endl;
    k63mer_to_graphposition.insert(kidxv.begin(),kidxv.end());
    sglib::OutputLog(sglib::INFO)<<"map ready"<<std::endl;
}
