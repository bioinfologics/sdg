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

    WorkSpace &ws;


};


#endif //BSG_UNTANGLER_HPP
