//
// Created by Bernardo Clavijo (EI) on 18/10/2017.
//

#include <iostream>
#include <fstream>
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
    std::string fasta_filename;
    //check the filename ends in .gfa
    if (filename.size()>4 and filename.substr(filename.size()-4,4)==".gfa"){
        fasta_filename=filename.substr(0,filename.size()-4)+".fasta";
    }
    else throw std::invalid_argument("filename of the gfa input does not end in gfa, it ends in '"+filename.substr(filename.size()-4,4)+"'");



    //load all sequences from fasta file if they're not canonical, flip and remember they're flipped
    std::cout<<"Loading sequences from "<<fasta_filename<<std::endl;
    std::ifstream fastaf(fasta_filename);
    if (!fastaf) throw std::invalid_argument("Can't read fasta file");
    std::string name,seq="",line;
    seq.reserve(10000000); //stupid hack but probably useful to reserve
    std::unordered_map<std::string,sgNodeID_t> oldnames_to_ids;
    nodes.clear();
    nodes.emplace_back(Node("")); //an empty node on 0, just to skip the space
    sgNodeID_t nextid=1;
    uint64_t rcnodes=0;
    while(!fastaf.eof()){
        std::getline(fastaf,line);
        if (fastaf.eof() or line[0]=='>'){
            if (!name.empty()) {
                //rough ansi C and C++ mix but it works
                nodes.emplace_back(Node(seq)); //an empty node on 0, just to skip the space
                if (oldnames_to_ids.find(name) != oldnames_to_ids.end())
                    throw std::logic_error("sequence " + name + " is already defined");
                oldnames_to_ids[name] = nextid++;
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
}

void SequenceGraph::write_to_gfa(std::string filename){
    std::string fasta_filename;
    //check the filename ends in .gfa
    if (filename.size()>4 and filename.substr(filename.size()-4,4)==".gfa"){
        fasta_filename=filename.substr(0,filename.size()-4)+".fasta";
    }
    else throw std::invalid_argument("filename of the gfa input does not end in gfa, it ends in '"+filename.substr(filename.size()-4,4)+"'");



    //load all sequences from fasta file if they're not canonical, flip and remember they're flipped
    std::cout<<"Loading sequences from "<<fasta_filename<<std::endl;
    std::ofstream fastaf(fasta_filename);
    if (!fastaf) throw std::invalid_argument("Can't write to fasta file");

    for (sgNodeID_t i=1;i<nodes.size();++i){
        fastaf<<">"<<i<<std::endl<<nodes[i].sequence<<std::endl;
    }
}