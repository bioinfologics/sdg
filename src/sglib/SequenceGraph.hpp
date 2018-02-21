//
// Created by Bernardo Clavijo (EI) on 18/10/2017.
//

#ifndef SG_SEQUENCEGRAPH_HPP
#define SG_SEQUENCEGRAPH_HPP

#include <vector>
#include <algorithm>
#include <string>
#include <map>
#include <unordered_map>
#include <set>
#include "sglib/readers/FileReader.h"

typedef int64_t sgNodeID_t; //first node is 1; negatives are RC

enum sgNodeStatus_t {sgNodeActive,sgNodeDeleted};

class Node{
public:
    Node(std::string _seq) : sequence(_seq),status(sgNodeActive){};
    std::string sequence;
    sgNodeStatus_t status;
    bool is_canonical();
    void make_rc();
};

class Link{
public:
    Link( sgNodeID_t _src, sgNodeID_t _dst, int32_t _dist) : source(_src), dest(_dst), dist(_dist) {};
    Link( const Link& l ) {
        source = l.source;
        dest = l.dest;
        dist = l.dist;
    };
    sgNodeID_t source,dest;
    int32_t dist;

    bool operator==(const  Link);
    bool operator<(const Link)const;

    friend std::ostream& operator<<(std::ostream& s, Link& l) {
        s << "Link from: " << l.source << " to: " << l.dest << " with a distance of: " << l.dist;
        return s;
    }

    Link make_canonical() const {
        Link rl(*this);
        std::swap(rl.source, rl.dest);
        return source < rl.source ? *this : rl;
    }

};

class SequenceGraphPath;

class SequenceGraph {
public:
    SequenceGraph(){};
    //=== I/O functions ===
    void load_from_gfa(std::string filename);
    void write_to_gfa(std::string filename);

    //=== graph operations ===
    sgNodeID_t add_node(Node n);
    void add_link( sgNodeID_t source, sgNodeID_t dest, int32_t d);

    std::vector<Link> get_fw_links( sgNodeID_t n);
    inline std::vector<Link> get_bw_links( sgNodeID_t n){ return get_fw_links (-n); };

    /*
     * Connected components, (TODO) optionally breaking up in repeats, nodes that class as repeats will be returned on their own
     */
    std::vector<std::vector<sgNodeID_t>> connected_components (int max_nr_totalinks=0, int max_nr_dirlinks=0, int min_rsize=0); //TODO: --> enable extra breaks in repeats

    // find bubbles in component of graph
    std::vector<std::vector<sgNodeID_t >> find_bubbles(std::vector<sgNodeID_t>);

    // remove_node
    void remove_node(sgNodeID_t);
    // remove_link
    void remove_link(sgNodeID_t source, sgNodeID_t dest);
    //These two need to mark expanded edges, and transfer read maps and unique kmers for non-expanded, but just read map for expanded.

    void join_path(SequenceGraphPath p, bool consume_nodes=true);
    // expand_path --> creates an edge with the consensus of a path, eliminates old nodes if only in path and unused edges
    void join_all_unitigs();
    std::vector<SequenceGraphPath> get_all_unitigs(uint16_t min_nodes);
    // simplify --> executes expand_path on every multi-sequence unitig


    // tip_clip -> eliminates tips.


    //void explode_node( sgNodeID_t node, uint16_t k);
    //void explode_all_nodes ();
    //void collapse_identical_nodes ();

    //later
    //project spectra and use for flow

    bool link_exists(sgNodeID_t from, sgNodeID_t to) {
        // Look for link between starting node and the new node.
        auto l = links[std::abs(from)].begin();
        // TODO: Can this just be a std::find?
        for (; l != links[std::abs(from)].end(); ++l){
            if (l->source == from and l->dest == to) break;
        }
        return l != links[std::abs(from)].end();
    }


    //=== internal variables ===
    std::vector<sgNodeID_t> oldnames_to_nodes(std::string _oldnames);
    std::vector<Node> nodes;
    std::vector<std::vector<Link>> links;
    std::unordered_map<std::string,sgNodeID_t> oldnames_to_ids;
    std::vector<std::string> oldnames;
    std::string filename,fasta_filename;
};


class SequenceGraphPath {
public:
    std::vector<sgNodeID_t> nodes;
    explicit SequenceGraphPath(SequenceGraph & _sg, std::vector<sgNodeID_t> _nodes={}) : sg(_sg) ,nodes(_nodes) {};
    std::string get_fasta_header() const;
    std::string get_sequence() const;
    void reverse();
    bool is_canonical();
    bool operator==(const SequenceGraphPath& rhs) const;
    bool operator<(const SequenceGraphPath& rhs) const;
    std::set<sgNodeID_t> make_set_of_nodes() const;

    bool append_to_path(sgNodeID_t newnode) {
        if (nodes.empty()) {
            //std::cout << "Path is new and empty." << std::endl;
            nodes.emplace_back(newnode);
            //std::cout << "Size of path is now: " << nodes.size() << std::endl;
            return true;
        }
        sgNodeID_t from = -(nodes.back());
        //std::cout << "Path is not new and empty: trying to append." << std::endl;
        bool append_possible = sg.link_exists(from, newnode);
        if (append_possible) {
            //std::cout << "Append is possible, adding node to the back." << std::endl;
            nodes.emplace_back(newnode);
            //std::cout << "Size of path is now: " << nodes.size() << std::endl;
        }
        return append_possible;
    }

    std::vector<Link> collect_links() const {
        std::vector<Link> s;
        sgNodeID_t pnode = 0;
        // just iterate over every node in path - contig names are converted to ids at construction
        for (auto &n:nodes) {
            if (pnode != 0) {
                // find link between pnode' output (+pnode) and n's sink (-n)
                auto l = sg.links[std::abs(pnode)].begin();
                for (; l != sg.links[std::abs(pnode)].end(); ++l){
                    if (l->source == pnode and l->dest == n) break;
                }
                if (l == sg.links[std::abs(pnode)].end()) {
                    std::cout << "can't find a link between " << pnode << " and " << -n << std::endl;
                    throw std::runtime_error("path has no link");
                } else {
                    s.emplace_back(*l);
                }
            }
            pnode=-n;
        }
        return s;
    }

private:
    SequenceGraph& sg;
};

class SequenceSubGraph {
public:
    std::vector<sgNodeID_t> nodes;
    explicit SequenceSubGraph(SequenceGraph & _sg, std::vector<sgNodeID_t> _nodes={})  : sg(_sg) ,nodes(_nodes) {};
    SequenceGraphPath make_path(); //returns empty path if not linear

private:
    SequenceGraph& sg;
};

#endif //SG_SEQUENCEGRAPH_HPP
