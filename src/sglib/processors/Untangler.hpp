//
// Created by Bernardo Clavijo (EI) on 26/02/2018.
//

#ifndef BSG_UNTANGLER_HPP
#define BSG_UNTANGLER_HPP


#include <sglib/WorkSpace.hpp>

class Untangler {
public:
    explicit Untangler(WorkSpace & _ws): ws(_ws) {};
    uint64_t solve_canonical_repeats_by_tags(std::unordered_set<uint64_t> & reads_to_remap);
    std::vector<std::pair<sgNodeID_t, sgNodeID_t>> get_all_HSPNPs();
    uint64_t extend_HSPNPs_by_tagwalking();

    WorkSpace &ws;


    std::vector<SequenceGraphPath> make_parallel_paths(std::vector<SequenceGraphPath>);

    bool all_nodes_consumed(std::vector<SequenceGraphPath>);

    std::vector<sgNodeID_t> shared_nodes(std::vector<std::vector<SequenceGraphPath>> parallel_paths);

    std::vector<SequenceGraphPath> combine( std::vector<SequenceGraphPath> parallel_paths1, std::vector<SequenceGraphPath> parallel_paths2 );
};


#endif //BSG_UNTANGLER_HPP
