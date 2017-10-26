//
// Created by Bernardo Clavijo (EI) on 18/10/2017.
//

#include <iostream>
#include <fstream>
#include <sstream>
#include "SequenceGraph.hpp"

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

void SequenceGraph::load_from_gfa(std::string filename) {
    std::string fasta_filename,line;
    //check the filename ends in .gfa
    if (filename.size()>4 and filename.substr(filename.size()-4,4)==".gfa"){
        fasta_filename=filename.substr(0,filename.size()-4)+".fasta";
    }
    else throw std::invalid_argument("filename of the gfa input does not end in gfa, it ends in '"+filename.substr(filename.size()-4,4)+"'");

    std::ifstream gfaf(filename);
    if (!gfaf) throw std::invalid_argument("Can't read gfa file");
    std::getline(gfaf, line);
    if (line!="H\tVN:Z:1.0") std::cout<<"WARNING, first line of gfa doesn't correspond to GFA1"<<std::endl;

    std::ifstream fastaf(fasta_filename);
    if (!fastaf) throw std::invalid_argument("Can't read fasta file");



    //load all sequences from fasta file if they're not canonical, flip and remember they're flipped
    std::cout<<"Loading sequences from "<<fasta_filename<<std::endl;

    std::string name,seq="";
    seq.reserve(10000000); //stupid hack but probably useful to reserve
    oldnames_to_ids.clear();
    nodes.clear();
    links.clear();
    add_node(Node("")); //an empty node on 0, just to skip the space
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

    std::string rtype,source,sourcedir,dest,destdir,cigar;
    sgNodeID_t src_id,dest_id;
    uint64_t lcount=0;
    Link l(0,0,0);
    while(std::getline(gfaf, line) and !gfaf.eof()) {
        //std::cout<<line<<std::endl;
        std::istringstream iss(line);
        iss >> rtype;
        if (rtype == "L"){
            iss >> source;
            iss >> sourcedir;
            iss >> dest;
            iss >> destdir;
            iss >> cigar;
            //std::cout<<"'"<<source<<"' '"<<sourcedir<<"' '"<<dest<<"' '"<<destdir<<"'"<<std::endl;
            if (oldnames_to_ids.find(source)==oldnames_to_ids.end()){
                oldnames_to_ids[source] = add_node(Node(""));
                //std::cout<<"added source!" <<source<<std::endl;
            }
            if (oldnames_to_ids.find(dest)==oldnames_to_ids.end()){
                oldnames_to_ids[dest] = add_node(Node(""));
                //std::cout<<"added dest! "<<dest<<std::endl;
            }
            src_id=oldnames_to_ids[source];
            dest_id=oldnames_to_ids[dest];

            //Go from GFA's "Links as paths" to a normal "nodes have sinks (-/start/left) and sources (+/end/right)
            if (src_id<0) {
                src_id = -src_id;
            }
            else {
                if (sourcedir=="-") sourcedir="+";
                else sourcedir="-";
            }

            if (dest_id<0) {
                dest_id=-dest_id;
                if (destdir=="-") destdir="+";
                else destdir="-";
            }

            if (sourcedir=="-") l.source=-src_id;
            else l.source=src_id;
            if (destdir=="-") l.dest=-dest_id;
            else l.dest=dest_id;
            l.dist=0;
            if (cigar.size()>1 and cigar[cigar.size()-1]=='M') {
                //TODO: better checks for M
                l.dist=-atoi(cigar.c_str());
                //std::cout<<l.dist<<std::endl;
            }
            links[src_id].emplace_back(l);
            std::swap(l.source,l.dest);
            links[dest_id].emplace_back(l);
            ++lcount;
        }
    }
    std::cout<<nodes.size()-1<<" nodes after connecting with "<<lcount<<" links"<<std::endl;
}

sgNodeID_t SequenceGraph::add_node(Node n) {
    nodes.emplace_back(n);
    links.emplace_back();
    return (sgNodeID_t) nodes.size()-1;
}

void SequenceGraph::write_to_gfa(std::string filename){
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
        fastaf<<">seq"<<i<<std::endl<<nodes[i].sequence<<std::endl;
        gfaf<<"S\tseq"<<i<<"\t*\tLN:i:"<<nodes[i].sequence.size()<<"\tUR:Z:"<<fasta_filename<<std::endl;
    }

    for (auto &ls:links){
        for (auto &l:ls)
            if (l.source<=l.dest) {
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
                std::cout<<"can't find a link between "<<pnode<<" and "<<-n<<std::endl;
                throw std::runtime_error("path has no link");
            } else {
                std::cout<<"validating distance between "<<pnode<<" and "<<n<<" of "<<l->dist<<std::endl;
                if (l->dist>0){
                    for (auto c=l->dist;c>0;--c) s+="N";
                }
                else {
                    auto ovl=-l->dist;
                    for (auto s1=s.c_str()+s.size()-ovl,s2=nseq.c_str();*s1!=NULL;++s1,++s2)
                        if (*s1!=*s2)
                            throw std::runtime_error("path overlap is invalid!");
                    std::cout<<"Overlap validated!"<<std::endl;
                    nseq.erase(0,ovl);
                }
            }
        }
        s+=nseq;
        pnode=-n;
    }
    return s;
}