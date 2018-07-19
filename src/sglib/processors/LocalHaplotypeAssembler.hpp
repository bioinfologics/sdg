//
// Created by Bernardo Clavijo (EI) on 18/07/2018.
//

#ifndef BSG_LOCALHAPLOTYPEASSEMBLER_HPP
#define BSG_LOCALHAPLOTYPEASSEMBLER_HPP


#include <sglib/WorkSpace.hpp>

class LocalHaplotypeAssembler {
public:
    //Constructor from ws and backbone, creates lists of tags, read ids, anchor sequences and whatnot.
    LocalHaplotypeAssembler(WorkSpace & _ws, std::vector<sgNodeID_t > _backbone);

    //TODO: Constructor from ws and problem file. Reads everything but the ws/graph from the file;
    LocalHaplotypeAssembler(WorkSpace & _ws, std::string problem_filename);

    //TODO: Constructor from full file. Reads everything including a slimmed down ws/graph from the file;
    LocalHaplotypeAssembler(std::string problem_filename);

    //TODO: write down to disk (just the problem, will need ws when loading)
    void write_problem(std::string prefix);

    //TODO: write down to disk (full set of things, including a slimmed down ws)
    void write_full(std::string prefix,bool keep_full_sg=false);

    //TODO: and yes, also perform the assembly ;)
    void path_all_reads(); //creates a path for every read, both 10x and LMP
    void assemble(int k, int min_cov, bool tag_cov);

    void patch_graph_in_workspace();

    WorkSpace & ws;
    std::vector<sgNodeID_t > backbone;
    std::vector<Node > backbone_nodes;
    std::set<bsg10xTag> tagSet;
    std::vector<std::pair<uint16_t , std::vector<uint64_t>>> paired_reads;


};


#endif //BSG_LOCALHAPLOTYPEASSEMBLER_HPP
