//
// Created by Luis Yanes (EI) on 22/03/2018.
//

#include <sdglib/graph/SequenceDistanceGraph.hpp>

std::string SequenceGraphPath::get_fasta_header(bool use_oldnames) const {
    std::string h = ">sgPath_[";
    if (use_oldnames) {
        for (auto &n:nodes) {
            h += n > 0 ? '+' : '-';
            h += sg.nodeID_to_name(std::abs(n)) + ", ";
        }
    } else {
        for (auto &n:nodes) {
            h += std::to_string(n)+", ";
        }
    }
    h.resize(h.size()-1);
    h += ']';
    return h;
}

std::string SequenceGraphPath::get_sequence() const {
    sdglib::OutputLog(sdglib::LogLevels::DEBUG) << "get_sequence for a SequenceGraphPath with nodes = [" ;
    if(sdglib::OutputLogLevel == sdglib::LogLevels::DEBUG) {
        for(auto &n:nodes) std::cout << " " << n;
        std::cout << " ]" << std::endl;
    }
    std::string s = "";
    sgNodeID_t pnode = 0;
    // just iterate over every node in path - contig names are converted to ids at construction
    for (auto &n:nodes) {
        std::string nseq;
        if (n>0){
            nseq = sg.nodes[n].sequence;
        } else {
            auto rcn = sg.nodes[-n];
            rcn.make_rc();
            nseq = rcn.sequence;
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
                if (l->dist > 0){
                    for (auto c = l->dist; c > 0; --c) s += "N";
                }
                else {
                    auto ovl =- l->dist;
                    for (auto s1 = s.c_str() + s.size() - ovl, s2 = nseq.c_str(); *s1 != NULL; ++s1, ++s2)
                        if (*s1 != *s2)
                            throw std::runtime_error("path overlap is invalid!");
                    nseq.erase(0, ovl);
                }
            }
        }
        s+=nseq;
        pnode=-n;
    }
    sdglib::OutputLog(sdglib::LogLevels::DEBUG) << "get_sequence finished successfully for SequenceGraphPath with nodes = [" ;
    if(sdglib::OutputLogLevel == sdglib::LogLevels::DEBUG) {
        for(auto &n:nodes) std::cout<<" "<<n;
        std::cout<<" ]"<<std::endl;
    }
    return s;
}

void SequenceGraphPath::reverse(){
    std::vector<sgNodeID_t> newn;
    for (auto n=nodes.rbegin();n<nodes.rend();++n) newn.emplace_back(-*n);
    nodes=newn;
}

bool SequenceGraphPath::is_canonical() {
    auto rp=*this;
    rp.reverse();
    return this->get_sequence()<rp.get_sequence();
}

std::set<sgNodeID_t> SequenceGraphPath::make_set_of_nodes() const {
    std::set<sgNodeID_t> s;
    std::transform(nodes.begin(),
                   nodes.end(),
                   std::inserter(s, s.end()),
                   [](const sgNodeID_t& n) -> sgNodeID_t { return std::abs(n); });
    return s;
}

bool SequenceGraphPath::operator==(const SequenceGraphPath& rhs) const {
    return make_set_of_nodes() == rhs.make_set_of_nodes();
}

bool SequenceGraphPath::operator<(const SequenceGraphPath& rhs) const {
    return make_set_of_nodes() < rhs.make_set_of_nodes();
}

SequenceGraphPath& SequenceGraphPath::operator=(const SequenceGraphPath &other) {
    if (&sg != &other.sg) { throw std::runtime_error("Can only copy paths from the same SequenceDistanceGraph"); }
    if (&other == this) {
        return *this;
    }
    nodes = other.nodes;
    return *this;
}

std::vector<Link> SequenceGraphPath::get_next_links() {
    return sg.get_fw_links(nodes.back());
}

size_t SequenceGraphPath::get_sequence_size_fast() const{
    size_t size=0;
    //std::string s="";
    sgNodeID_t pnode=0;
    // just iterate over every node in path - contig names are converted to ids at construction
    for (auto &n:nodes) {
        std::string nseq;
        size+=sg.nodes[llabs(n)].sequence.size();
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
