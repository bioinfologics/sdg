//
// Created by Bernardo Clavijo (EI) on 18/10/2017.
//


#include "SequenceDistanceGraph.hpp"
#include <iostream>
#include <fstream>
#include <sstream>
#include <set>
#include <cmath>
#include <list>
#include <queue>
#include <stack>
#include <tuple>
#include <functional>
#include <sdglib/utilities/io_helpers.hpp>

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
                sdglib::OutputLog(sdglib::LogLevels::WARN) << "Unexpected character in fasta file: '" << r << "'" << std::endl;
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

std::string SequenceDistanceGraph::get_node_sequence(sgNodeID_t n){
    if (n>0) return nodes[n].sequence;
    n=-n;
    std::string rseq;
    rseq.resize(nodes[n].sequence.size());
    for (size_t i=0,j=nodes[n].sequence.size()-1;i<nodes[n].sequence.size();++i,--j){
        switch(nodes[n].sequence[j]){
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
    return rseq;
}

uint64_t SequenceDistanceGraph::get_node_size(sgNodeID_t n) {
    return nodes[llabs(n)].sequence.size();
}

bool SequenceDistanceGraph::is_sane() const {
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
    return true;
}

sgNodeID_t SequenceDistanceGraph::add_node(Node n) {
    nodes.emplace_back(n);
    links.emplace_back();
    return (sgNodeID_t) nodes.size()-1;
}

void SequenceDistanceGraph::remove_node(sgNodeID_t n) {
    sgNodeID_t node=(n>0? n:-n);
    auto oldlinks=links[node];//this creates a copy to allow the iteration
    for (auto &l:oldlinks) remove_link(l.source,l.dest);
    nodes[node].status=NodeStatus::Deleted;
    //TODO: this is a lazy solution
    nodes[node].sequence.clear();
    //TODO: remove read mappings
}



size_t SequenceDistanceGraph::count_active_nodes() {
    size_t t = 0;
    for (auto &n: nodes) {
        if (n.status == NodeStatus::Active) ++t;
    }
    return t;
}

bool Link::operator==(const Link a) const {
    if (a.source == source && a.dest == dest){
        return true;
    }
    return false;
}

bool Link::operator<(const Link a) const {
    if (a.source < source){
        return true;
    } if (a.source == source && a.dest < dest) {
        return  true;
    }
    return false;
}


std::vector<std::vector<sgNodeID_t>> SequenceDistanceGraph::connected_components(int max_nr_totalinks, int max_nr_dirlinks,
                                                                         int min_rsize) {
    std::vector<bool> used(nodes.size());
    std::vector<std::vector<sgNodeID_t>> components;
    //TODO: first find all repeats, add them as independent components and mark them as used.
    size_t max_component = 0;
    for (sgNodeID_t start_node=1;start_node<nodes.size();++start_node){
        if (nodes[start_node].status==NodeStatus::Deleted) continue;
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

void SequenceDistanceGraph::write(std::ofstream & output_file) {
    uint64_t count;
    count=nodes.size();

    output_file.write((char *) &count,sizeof(count));
    for (auto &n:nodes){
        output_file.write((char *) &n.status,sizeof(n.status));
        count=n.sequence.size();
        output_file.write((char *) &count,sizeof(count));
        output_file.write((char *) n.sequence.data(),count);
    }

    sdglib::write_flat_vectorvector(output_file, links);
}

void SequenceDistanceGraph::read(std::ifstream & input_file) {
    uint64_t count;
    input_file.read((char *) &count,sizeof(count));
    nodes.clear();
    nodes.reserve(count);
    uint64_t active=0;
    for (auto i=0;i<count;++i){
        uint64_t seqsize;
        std::string seq;
        NodeStatus status;
        input_file.read((char *) &status,sizeof(status));
        input_file.read((char *) &seqsize,sizeof(seqsize));
        seq.resize(seqsize);
        input_file.read((char *) seq.data(),seqsize);
        nodes.emplace_back(seq);
        nodes.back().status=status;
        if (nodes.back().status==NodeStatus::Active) ++active;
    }

    sdglib::read_flat_vectorvector(input_file, links);
}

std::string SequenceDistanceGraph::ls(int level, bool recursive) {
    std::stringstream ss;
    std::string spacer(2 * level, ' ');
    uint64_t linkcount = 0;
    for (auto lv:links)
        for (auto l:lv)
            if (llabs(l.source) <= llabs(l.dest)) ++linkcount;
    ss << spacer << "SequenceDistanceGraph " << name <<": "<< count_active_nodes() <<" nodes, " << linkcount << " links" << std::endl;
    return ss.str();
}

void SequenceDistanceGraph::load_from_gfa(std::string filename) {
    std::string line;
    filename=filename;
    //check the filename ends in .gfa
    if (filename.size()>4 and filename.substr(filename.size()-4,4)==".gfa"){
        fasta_filename=filename.substr(0,filename.size()-4)+".fasta";
    }
    else if (filename.size()>5 and filename.substr(filename.size()-5,5)==".gfa2"){
        fasta_filename=filename.substr(0,filename.size()-5)+".fasta";
    }
    else throw std::invalid_argument("Filename of the gfa input does not end in gfa, it ends in '"+filename.substr(filename.size()-4,4)+"'");

    std::ifstream gfaf(filename);
    if (!gfaf) {
        std::cerr << "Failed to open " << filename <<": " << strerror(errno);
        throw std::invalid_argument("Can't read gfa file");
    }
    if (gfaf.peek() == std::ifstream::traits_type::eof()) throw std::invalid_argument("Empty gfa file");

    std::ifstream fastaf(fasta_filename);
    if (!fastaf) {
        fasta_filename=fasta_filename.substr(0,fasta_filename.size()-6)+".fa";
        fastaf=std::ifstream(fasta_filename);
    }
    if (!fastaf) {
        std::cerr << "Failed to open " << fasta_filename <<": " << strerror(errno);
        throw std::invalid_argument("Can't read graph fasta file");
    }
    if (fastaf.peek() == std::ifstream::traits_type::eof()) throw std::invalid_argument("Empty fasta file");

    std::getline(gfaf, line);
    if (line=="H\tVN:Z:1.0") {
        load_from_gfa1(gfaf, fastaf);
    } else if (line == "H\tVN:Z:2.0") {
        load_from_gfa2(gfaf, fastaf);
    } else {
        std::cout<<"WARNING, unsuported GFA version defined on header: " << line << std::endl;
        std::cout<<"The file will be treated as GFA1" << std::endl;
        load_from_gfa1(gfaf, fastaf);
    }
}

void SequenceDistanceGraph::load_from_fasta(std::string filename) {
    std::string line;
    filename=filename;
    fasta_filename=filename;

    std::ifstream fastaf(fasta_filename);
    sdglib::OutputLog(sdglib::LogLevels::INFO) << "Graph fasta filesname: " << fasta_filename << std::endl;
    if (!fastaf) throw std::invalid_argument("Can't read graph fasta file");
    if (fastaf.peek() == std::ifstream::traits_type::eof()) throw std::invalid_argument("Empty fasta file");

    //load all sequences from fasta file if they're not canonical, flip and remember they're flipped
    sdglib::OutputLog(sdglib::LogLevels::INFO) << "Loading sequences from " << fasta_filename << std::endl;

    std::string name, seq = "";
    seq.reserve(10000000); //stupid hack but probably useful to reserve
    oldnames_to_ids.clear();
    oldnames.push_back("");
    nodes.clear();
    links.clear();
    add_node(Node("",NodeStatus::Deleted)); //an empty deleted node on 0, just to skip the space
    uint64_t rcnodes=0;
    while(!fastaf.eof()){
        std::getline(fastaf,line);
        if (fastaf.eof() or line[0] == '>'){

            if (!name.empty()) {
                //rough ansi C and C++ mix but it works
                if (oldnames_to_ids.find(name) != oldnames_to_ids.end())
                    throw std::logic_error("sequence " + name + " is already defined");
                oldnames_to_ids[name] = add_node(Node(seq));
                oldnames.push_back(name);
            }

            // Clear the name and set name to the new name, this is a new sequence!
            name.clear();
            for (auto i = 1; i < line.size() and line[i] != ' '; ++i) name += line[i];
            seq = "";
        } else {
            seq += line;

        }
    }
    sdglib::OutputLog(sdglib::LogLevels::INFO) << nodes.size()-1 << " nodes loaded " << std::endl;
}

std::vector<sgNodeID_t> SequenceDistanceGraph::oldnames_to_nodes(std::string _oldnames) {
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

std::vector<SequenceGraphPath> SequenceDistanceGraph::get_all_unitigs(uint16_t min_nodes) {
    std::vector<SequenceGraphPath> unitigs;
    std::vector<bool> used(nodes.size(),false);

    for (auto n=1;n<nodes.size();++n){
        if (used[n] or nodes[n].status==NodeStatus::Deleted) continue;
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

uint32_t SequenceDistanceGraph::join_all_unitigs() {
    uint32_t joined=0;
    for (auto p:get_all_unitigs(2)){
        join_path(p);
        ++joined;
    }
    return joined;
}

void SequenceDistanceGraph::join_path(SequenceGraphPath p, bool consume_nodes) {
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


void SequenceDistanceGraph::expand_node(sgNodeID_t nodeID, std::vector<std::vector<sgNodeID_t>> bw,
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

std::vector<SequenceSubGraph> SequenceDistanceGraph::get_all_bubbly_subgraphs(uint32_t maxsubgraphs) {
    std::vector<SequenceSubGraph> subgraphs;
    std::vector<bool> used(nodes.size(),false);
    /*
     * the loop always keep the first and the last elements as c=2 collapsed nodes.
     * it starts with a c=2 node, and goes thorugh all bubbles fw, then reverts the subgraph and repeats
     */
    SequenceSubGraph subgraph(*this);
    for (auto n=1;n<nodes.size();++n){
        if (used[n] or nodes[n].status==NodeStatus::Deleted) continue;
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
            //for (auto &n:subgraph) std::cout<<"  "<<n<<" ("<<sdg.nodes[(n>0?n:-n)].sequence.size()<<"bp)";
            //std::cout<<std::endl;
        }
    }
    return subgraphs;
}

void SequenceDistanceGraph::print_bubbly_subgraph_stats(const std::vector<SequenceSubGraph> &bubbly_paths) {
    std::vector<uint64_t> solved_sizes,original_sizes;
    sdglib::OutputLog()<<bubbly_paths.size()<<" bubbly paths"<<std::endl;
    auto &log_no_date=sdglib::OutputLog(sdglib::LogLevels::INFO,false);
    uint64_t total_size=0,total_solved_size=0;
    for (const auto &bp:bubbly_paths){
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
    sdglib::OutputLog() <<"Currently "<<original_sizes.size()<<" sequences with "<<total_size<<"bp, ";
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
    sdglib::OutputLog()<<"Potentially "<<solved_sizes.size()<<" sequences with "<<total_solved_size<<"bp, ";
    acc=0;
    for (auto s:solved_sizes){
        if (acc<sn20s and acc+s>sn20s) log_no_date<<"N20: "<<s<<"  ";
        if (acc<sn50s and acc+s>sn50s) log_no_date<<"N50: "<<s<<"  ";
        if (acc<sn80s and acc+s>sn80s) log_no_date<<"N80: "<<s<<"  ";
        acc+=s;
    }
    log_no_date<<std::endl;

}

void SequenceDistanceGraph::print_status() {
    auto &log_no_date=sdglib::OutputLog(sdglib::LogLevels::INFO,false);
    std::vector<uint64_t> node_sizes;
    uint64_t total_size=0;
    for (auto n:nodes) {
        if (n.status!=NodeStatus::Deleted) {
            total_size += n.sequence.size();
            node_sizes.push_back(n.sequence.size());
        }
    }
    std::sort(node_sizes.rbegin(),node_sizes.rend());
    auto on20s=total_size*.2;
    auto on50s=total_size*.5;
    auto on80s=total_size*.8;
    sdglib::OutputLog() <<"The graph contains "<<node_sizes.size()<<" sequences with "<<total_size<<"bp, ";
    uint64_t acc=0;
    for (auto s:node_sizes){
        if (acc==0)  log_no_date<<"N0: "<<s<<"bp  ";
        if (acc<on20s and acc+s>on20s) log_no_date<<"N20: "<<s<<"bp  ";
        if (acc<on50s and acc+s>on50s) log_no_date<<"N50: "<<s<<"bp  ";
        if (acc<on80s and acc+s>on80s) log_no_date<<"N80: "<<s<<"bp  ";
        acc+=s;
        if (acc==total_size)  log_no_date<<"N100: "<<s<<"bp  ";
    }
    log_no_date<<std::endl;
}

void SequenceDistanceGraph::load_from_gfa1(std::ifstream &gfaf, std::ifstream &fastaf) {
    std::string line;
    sdglib::OutputLog(sdglib::LogLevels::INFO) << "Graph fasta filesname: " << fasta_filename << std::endl;
    //load all sequences from fasta file if they're not canonical, flip and remember they're flipped
    sdglib::OutputLog(sdglib::LogLevels::INFO) << "Loading sequences from " << fasta_filename << std::endl;

    std::unordered_map<std::string, unsigned int> gap_dist;
    std::string name, seq = "";
    seq.reserve(10000000); //stupid hack but probably useful to reserve
    oldnames_to_ids.clear();
    oldnames.push_back("");
    nodes.clear();
    links.clear();
    add_node(Node("",NodeStatus::Deleted)); //an empty deleted node on 0, just to skip the space
    sgNodeID_t nextid=1;
    uint64_t rcnodes=0;
    while(!fastaf.eof()){
        std::getline(fastaf,line);
        if (fastaf.eof() or line[0] == '>'){

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

            // Clear the name and set name to the new name, this is a new sequence!
            name.clear();
            for (auto i = 1; i < line.size() and line[i] != ' '; ++i) name += line[i];
            seq = "";
        } else {
            seq += line;

        }
    }
    sdglib::OutputLog(sdglib::LogLevels::INFO) << nodes.size()-1 << " nodes loaded (" << rcnodes << " canonised)." << std::endl;

    //load store all connections.

    std::string gfa_rtype,gfa_source,gfa_sourcedir,gfa_dest,gfa_destdir,gfa_cigar,gfa_star,gfa_length,gap_id,gap_dir,seq_location;
    sgNodeID_t src_id,dest_id;
    int32_t dist;
    uint64_t lcount=0;
    uint64_t dist_egt0(0);
    while(std::getline(gfaf, line) and !gfaf.eof()) {
        std::istringstream iss(line);
        iss >> gfa_rtype;

        if (gfa_rtype == "S") {
            iss >> gfa_source;
            iss >> gfa_star;
            iss >> gfa_length; // parse to number
            if (!iss.eof()) iss >> seq_location;
            if (gfa_star != "*") {
                throw std::logic_error("Sequences should be in a separate file.");
            }

            if (gfa_star == "*" and seq_location.empty()) {
                gap_dist[gfa_source] = std::stoi(gfa_length.substr(5));
            } else {
                // Check equal length seq and node length reported in gfa
                if (oldnames_to_ids.find(gfa_source) != oldnames_to_ids.end()) {
                    if (std::stoi(gfa_length.substr(5)) !=
                        nodes[std::abs(oldnames_to_ids[gfa_source])].sequence.length()) {
                        throw std::logic_error(
                                "Different length in node and fasta for sequence: " + gfa_source + " -> gfa:" +
                                gfa_length.substr(5) + ", fasta: " +
                                std::to_string(nodes[oldnames_to_ids[gfa_source]].sequence.length()));
                    }
                }
            }
            seq_location.clear();
        } else if (gfa_rtype == "L"){
            iss >> gfa_source;
            iss >> gfa_sourcedir;
            iss >> gfa_dest;
            iss >> gfa_destdir;
            iss >> gfa_cigar;

            if (gap_dist.find(gfa_dest) != gap_dist.end()) {
                std::getline(gfaf, line);
                std::istringstream gap_ss(line);
                gap_ss >> gfa_rtype;
                gap_ss >> gap_id;        // Ignore gap_id x2
                gap_ss >> gap_dir;     // Ignore gap_dir x2
                gap_ss >> gfa_dest;        // Read dest
                gap_ss >> gfa_destdir;     // Read dest_dir
                gap_ss >> gfa_cigar;       // Ignore cigar
            }
            //std::cout<<"'"<<source<<"' '"<<gfa_sourcedir<<"' '"<<dest<<"' '"<<destdir<<"'"<<std::endl;
            if (gap_dist.find(gap_id) == gap_dist.end()) {
                if (oldnames_to_ids.find(gfa_source) == oldnames_to_ids.end()) {
                    oldnames_to_ids[gfa_source] = add_node(Node(""));
                    //std::cout<<"added source!" <<source<<std::endl;
                }
                if (oldnames_to_ids.find(gfa_dest) == oldnames_to_ids.end()) {
                    oldnames_to_ids[gfa_dest] = add_node(Node(""));
                    //std::cout<<"added dest! "<<dest<<std::endl;
                }
            }
            src_id=oldnames_to_ids[gfa_source];
            dest_id=oldnames_to_ids[gfa_dest];

            //Go from GFA's "Links as paths" to a normal "nodes have sinks (-/start/left) and sources (+/end/right)
            if (gfa_sourcedir == "+") src_id=-src_id;
            if (gfa_destdir == "-") dest_id=-dest_id;
            dist=0;
            if (gfa_cigar.size()>1 and gfa_cigar[gfa_cigar.size()-1]=='M') {
                //TODO: better checks for M
                dist=-atoi(gfa_cigar.c_str());
                if (dist>=0) {
                    dist_egt0++;
                }
            }

            if (!gap_id.empty()) {
                const auto d = gap_dist.find(gap_id);
                if (d != gap_dist.cend()) {
                    dist = d->second;
                } else {
                    throw std::runtime_error(gfa_dest + " has not been seen in " + filename +
                                             " yet, please ensure gaps are defined before being referenced");
                }
                gap_id.clear();
            }

            add_link(src_id,dest_id,dist);
            ++lcount;
        }
    }
    if (dist_egt0 > lcount*0.5f) {
        sdglib::OutputLog(sdglib::LogLevels::WARN) << "Warning: The loaded graph contains " << dist_egt0 << " non-overlapping links out of " << lcount << std::endl;
    }
    sdglib::OutputLog(sdglib::LogLevels::INFO) << nodes.size() - 1 << " nodes after connecting with " << lcount << " links." << std::endl;
}

void SequenceDistanceGraph::load_from_gfa2(std::ifstream &gfaf, std::ifstream &fastaf) {
    std::string line;
    sdglib::OutputLog(sdglib::LogLevels::INFO) << "Graph fasta filename: " << fasta_filename << std::endl;
    //load all sequences from fasta file if they're not canonical, flip and remember they're flipped
    sdglib::OutputLog(sdglib::LogLevels::INFO) << "Loading sequences from " << fasta_filename << std::endl;

    std::unordered_map<std::string, unsigned int> gap_dist;
    std::string name, seq = "";
    seq.reserve(10000000); //stupid hack but probably useful to reserve
    oldnames_to_ids.clear();
    oldnames.push_back("");
    nodes.clear();
    links.clear();
    add_node(Node("",NodeStatus::Deleted)); //an empty deleted node on 0, just to skip the space
    sgNodeID_t nextid=1;
    uint64_t rcnodes=0;
    while(!fastaf.eof()){
        std::getline(fastaf,line);
        if (fastaf.eof() or line[0] == '>'){

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

            // Clear the name and set name to the new name, this is a new sequence!
            name.clear();
            for (auto i = 1; i < line.size() and line[i] != ' '; ++i) name += line[i];
            seq = "";
        } else {
            seq += line;

        }
    }
    sdglib::OutputLog(sdglib::LogLevels::INFO) << nodes.size()-1 << " nodes loaded (" << rcnodes << " canonised)." << std::endl;

    //load store all connections.

    std::string gfa_rtype,gfa_source,gfa_sourcedir,gfa_dest,gfa_destdir,gfa_cigar,gfa_star,gfa_length,id;
    unsigned int source_beg,source_end,dest_beg,dest_end;
    sgNodeID_t src_id,dest_id;
    int32_t dist;
    uint64_t lcount=0;
    uint64_t dist_egt0(0);
    while(std::getline(gfaf, line) and !gfaf.eof()) {
        std::istringstream iss(line);
        iss >> gfa_rtype;

        if (gfa_rtype == "S") {
            iss >> gfa_source;
            iss >> gfa_length;
            iss >> gfa_star;

            // The length of a sequence reported by GFA2 isn't required to match the actual length of the sequence
            /*
            if (oldnames_to_ids.find(gfa_source) != oldnames_to_ids.end()) {
                if (std::stoi(gfa_length.substr(5)) !=
                    nodes[std::abs(oldnames_to_ids[gfa_source])].sequence.length()) {
                    throw std::logic_error(
                            "Different length in node and fasta for sequence: " + gfa_source + " -> gfa:" +
                            gfa_length.substr(5) + ", fasta: " +
                            std::to_string(nodes[oldnames_to_ids[gfa_source]].sequence.length()));
                }
            }*/

        } else if (gfa_rtype == "E") {
            iss >> id;
            iss >> gfa_source;
            iss >> gfa_dest;
            iss >> source_beg;
            iss >> source_end;
            iss >> dest_beg;
            iss >> dest_end;
            iss >> gfa_cigar;

            dist =  -1 * (source_end - source_beg);
            gfa_sourcedir = gfa_source.back();
            gfa_source.pop_back();
            gfa_destdir = gfa_dest.back();
            gfa_dest.pop_back();

            //std::cout<<"'"<<source<<"' '"<<gfa_sourcedir<<"' '"<<dest<<"' '"<<destdir<<"'"<<std::endl;
            if (oldnames_to_ids.find(gfa_source) == oldnames_to_ids.end()) {
                oldnames_to_ids[gfa_source] = add_node(Node(""));
                //std::cout<<"added source!" <<source<<std::endl;
            }
            if (oldnames_to_ids.find(gfa_dest) == oldnames_to_ids.end()) {
                oldnames_to_ids[gfa_dest] = add_node(Node(""));
                //std::cout<<"added dest! "<<dest<<std::endl;
            }

            src_id=oldnames_to_ids[gfa_source];
            dest_id=oldnames_to_ids[gfa_dest];

            //Go from GFA's "Links as paths" to a normal "nodes have sinks (-/start/left) and sources (+/end/right)
            if (gfa_sourcedir == "+") src_id=-src_id;
            if (gfa_destdir == "-") dest_id=-dest_id;

            add_link(src_id,dest_id,dist);
            ++lcount;
        } else if (gfa_rtype == "G") {
            iss >> id;
            iss >> gfa_source;
            iss >> gfa_dest;
            iss >> dist;

            gfa_sourcedir = gfa_source.back();
            gfa_source.pop_back();
            gfa_destdir = gfa_dest.back();
            gfa_dest.pop_back();

            //std::cout<<"'"<<source<<"' '"<<gfa_sourcedir<<"' '"<<dest<<"' '"<<destdir<<"'"<<std::endl;
            if (oldnames_to_ids.find(gfa_source) == oldnames_to_ids.end()) {
                oldnames_to_ids[gfa_source] = add_node(Node(""));
                //std::cout<<"added source!" <<source<<std::endl;
            }
            if (oldnames_to_ids.find(gfa_dest) == oldnames_to_ids.end()) {
                oldnames_to_ids[gfa_dest] = add_node(Node(""));
                //std::cout<<"added dest! "<<dest<<std::endl;
            }

            src_id=oldnames_to_ids[gfa_source];
            dest_id=oldnames_to_ids[gfa_dest];

            //Go from GFA's "Links as paths" to a normal "nodes have sinks (-/start/left) and sources (+/end/right)
            if (gfa_sourcedir == "+") src_id=-src_id;
            if (gfa_destdir == "-") dest_id=-dest_id;

            add_link(src_id,dest_id,dist);

            dist_egt0++;
        }
    }
    if (dist_egt0 > lcount*0.5f) {
        sdglib::OutputLog(sdglib::LogLevels::WARN) << "Warning: The loaded graph contains " << dist_egt0 << " non-overlapping links out of " << lcount << std::endl;
    }
    sdglib::OutputLog(sdglib::LogLevels::INFO) << nodes.size() - 1 << " nodes after connecting with " << lcount << " links." << std::endl;

}
