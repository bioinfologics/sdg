//
// Created by Luis Yanes (EI) on 22/03/2018.
//
#include <string>
#include <fstream>
#include <sglib/graph/SequenceGraph.hpp>


void SequenceSubGraph::write_to_gfa(std::string filename) {
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
    std::unordered_set<sgNodeID_t > nodes_in_links;
    for (const auto &n:nodes) {
        for (auto &l:sg.links[std::abs(n)]) {
            nodes_in_links.insert(std::abs(l.source));
            nodes_in_links.insert(std::abs(l.dest));
        }
    }

    for (const auto &n:nodes_in_links){
        fastaf<<">"<<sg.oldnames[n]<<std::endl<<sg.nodes[std::abs(n)].sequence<<std::endl;
        gfaf<<"S\t"<<sg.oldnames[n]<<"\t*\tLN:i:"<<sg.nodes[std::abs(n)].sequence.length()<<"\tUR:Z:"<<fasta_filename<<std::endl;
    }

    for (const auto &n:nodes) {
        for (auto &l:sg.links[std::abs(n)]) {
            gfaf << "L\t";
            if (l.source > 0) gfaf << "" << sg.oldnames[std::abs(l.source)] << "\t-\t";
            else gfaf << "" << sg.oldnames[std::abs(l.source)] << "\t+\t";
            if (l.dest > 0) gfaf << "" << sg.oldnames[std::abs(l.dest)] << "\t+\t";
            else gfaf << "" << sg.oldnames[std::abs(l.dest)] << "\t-\t";
            gfaf << (l.dist < 0 ? -l.dist : 0) << "M" << std::endl;
        }
    }
}

uint64_t SequenceSubGraph::total_size() const {
    uint64_t t=0;
    for (auto &n:nodes) t+=sg.nodes[llabs(n)].sequence.size();
    return t;
}
