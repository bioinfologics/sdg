//
// Created by Bernardo Clavijo (EI) on 30/04/2020.
//

#pragma once
#include <sdglib/graph/DistanceGraph.hpp>
#include <sdglib/indexers/NKmerIndex.hpp>
#include <memory>
#include "LongReadsRecruiter.hpp"

/** @brief Representation of a part of a match in the PerfectMatchExtender (internal use)
 * it's a temporal storage for the match parts while the algorithm is constructing the paths across the overlaps of
 * the graph.
 *
 * The extend function here only extends within a node.
 */
class PerfectMatchPart{
public:
    /** @brief Extend a match for the read in the current node. does not jump nodes.
     *
     * Uses string matching to determine the match.
     *
     * @param readseq Read sequence
     * @param nodeseq Node sequence
     */
    void extend(const std::string & readseq,const std::string & nodeseq);


    sgNodeID_t node;

    uint64_t offset;//will only be set if first part
    int64_t previous_part;//will only be set if not first part
    uint64_t read_position;//position of the last matched base in read
    uint64_t node_position;//postiion of the last matched base in node (canonical orientation)

    bool completed_node=false;
    bool completed_read=false;
    bool extended=false;
    bool invalid=false;
};

/**
 * Extends a match starting from a kmer by following the matches even across links in successive nodes in the graph
 */
class PerfectMatchExtender{
public:
    PerfectMatchExtender(DistanceGraph & _dg, uint8_t _k):dg(_dg),k(_k){};

    /** @brief setup the read and reset structures
     *
     * @param _readseq read sequence to match and extend
     */
    void set_read(const std::string & _readseq);

    /** @brief reset the internal state of the class
     *
     */
    void reset();

    /** @brief setup the stating node to expand a match from
     *
     * @param node_id Node id
     * @param read_offset offset start in the read
     * @param node_offset offset start in the node
     */
    void add_starting_match(sgNodeID_t node_id, uint64_t read_offset, uint64_t node_offset);

    /** @brief Extend all possible matches FW jumping overlaping nodes.
     *
     * Matchparts are stored in the matchparts vector ( std::vector<PerfectMatchPart> )
     *
     * This function performs the extension following the topology of the provided graph, if a match is valid bridging
     * the connection between 2 nodes the match is cointinued and placed in the matchparts vector.
     * If there is no topology (dg with no links) this can be used but the extension will finish at the end of the node.
     */
    void extend_fw();

    /** @brief Sets the best_path to the best path of the matchparts collection
     *
     * Fills the std::vector<sgNodeID_t> best_path vector with a list of nodes in the path
     * The best path is defined as the path that goes further in the read.
     *
     * @param fill_offsets
     */
    void set_best_path(bool fill_offsets=false); //Todo: return pointer to the last part?

    /** @brief reconstructs the best path as a vector of perfect matches
     *
     * Fills std::vector<PerfectMatch> best_path_matches.
     *
     */
    void make_path_as_perfect_matches();


    DistanceGraph & dg;
    uint8_t k;
    uint32_t start_mp_readpos;
    std::vector<PerfectMatchPart> matchparts;
    std::vector<PerfectMatch> best_path_matches;
    std::vector<sgNodeID_t> best_path;
    std::vector<std::pair<uint32_t,uint32_t>> best_path_offsets;
    uint32_t best_path_offset;
    uint64_t last_readpos;
    uint64_t last_nodepos;
    std::string readseq;
    int64_t winning_last_part=-1;
    std::vector<int> votes;//optimisation, since best_path takes ages

};

/*class PerfectMatcher {
public:
    //Create with a graph and parameters for the index.
    PerfectMatcher(DistanceGraph &_dg,std::shared_ptr<NKmerIndex> _nki):dg(_dg),nki(_nki){};

    //Create with a graph and a pointer to the index.
    PerfectMatcher(DistanceGraph &_dg,uint8_t _k, uint16_t _max_freq):dg(_dg),nki(std::make_shared<NKmerIndex>(dg.sdg,_k,_max_freq)){};

    //set a sequence to map
    void set_sequence();

    //get_next_match returns a PerfectMatchPath object that points to the pme and enables to get the match parts
    PerfectMatchPart & get_next_part();

    void reset_part_getter()

    std::shared_ptr<NKmerIndex> nki;
    DistanceGraph & dg;
};*/