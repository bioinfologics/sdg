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
    uint64_t expand_canonical_repeats_by_tags(float min_ci, float max_ci, int min_tags=10);
    std::vector<std::pair<sgNodeID_t, sgNodeID_t>> get_all_HSPNPs();
    uint64_t extend_HSPNPs_by_tagwalking();

    void analise_paths_through_nodes();
    std::vector<std::pair<sgNodeID_t,sgNodeID_t>> find_bubbles(uint32_t min_size,uint32_t max_size);
    std::vector<std::pair<sgNodeID_t,sgNodeID_t>> solve_bubbly_paths();
    void pop_errors_by_ci_and_paths();

    std::vector<std::vector<std::pair<sgNodeID_t,uint32_t>>> find_tag_neighbours(uint32_t min_size, float min_ci, float max_ci);

    WorkSpace &ws;


    std::vector<SequenceGraphPath> make_parallel_paths(std::vector<SequenceGraphPath>);

    bool all_nodes_consumed(std::vector<SequenceGraphPath>);

    std::vector<sgNodeID_t> shared_nodes(std::vector<std::vector<SequenceGraphPath>> parallel_paths);

    std::vector<SequenceGraphPath> combine( std::vector<SequenceGraphPath> parallel_paths1, std::vector<SequenceGraphPath> parallel_paths2 );

    void connect_neighbours(uint64_t min_size, float min_ci, float max_ci, int64_t max_distance);
    void connect_neighbours_trivial(uint64_t min_size, float min_ci, float max_ci, int64_t max_distance,
                                     const std::vector<std::vector<std::pair<sgNodeID_t,uint32_t>>> &tagneighbours,
                                     const std::vector<std::vector<std::pair<sgNodeID_t,int64_t>>> & bndist,
                                     const std::vector<std::vector<std::pair<sgNodeID_t,int64_t>>> & fndist);

    void connect_neighbours_paths_to_same(uint64_t min_size, float min_ci, float max_ci, int64_t max_distance,
                                          const std::vector<std::vector<std::pair<sgNodeID_t,uint32_t>>> &tagneighbours,
                                          const std::vector<std::vector<std::pair<sgNodeID_t,int64_t>>> & bndist,
                                          const std::vector<std::vector<std::pair<sgNodeID_t,int64_t>>> & fndist);
};


#endif //BSG_UNTANGLER_HPP
