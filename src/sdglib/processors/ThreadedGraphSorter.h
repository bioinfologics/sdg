//
// Created by Gonzalo Garcia (EI) on 2020-11-10.
//

#ifndef SDG_THREADEDGRAPHSORTER_H
#define SDG_THREADEDGRAPHSORTER_H


#include <sdglib/views/NodeView.hpp>
#include <sdglib/mappers/LongReadsRecruiter.hpp>

class LocalOrder{
public:
    LocalOrder(){};
    explicit LocalOrder(const std::vector<sgNodeID_t> & nodes);
    std::vector<sgNodeID_t> as_signed_nodes() const;
    int64_t get_node_position(sgNodeID_t nid) const;
    LocalOrder reverse() const;
    void validate_with_rtg(const ReadThreadsGraph & rtg);
    void cleanup_with_rtg(const ReadThreadsGraph & rtg);
    auto size() const {return node_positions.size();};

    //explicit LocalOrder(std::vector<sgNodeID_t>);//construct a local order from a list of nodes
    LocalOrder merge(const LocalOrder & other,int max_overhang=4,float min_shared_perc=.5,int min_shared=20,float max_disordered_perc=.02) const;
    std::unordered_map<sgNodeID_t , int64_t > node_positions;
};

/** @brief Checks node happiness in an order
 * Compares the position in the proposed order with the links in the thread graph (trg_nt).
 *
 * The position of the node in the ordered is compared to the position of the node in each read thread, if the position
 * agrees the the node happines increases, if not decreases.
 *
 * The reads are represented as paths in the threads graph.
 *
 * Return an array of 3 scores [happy count, unhappy count, disconnected count]
 *  - Happy count is the number of reads where the order in the proposed order matches the order in the reads (threaded graph) both fw and bw.
 *  - The unhappy count is the number of reads were the order proposed don't match the order in the threading reads either fw or bw.
 *  - The disconnected count is the number of reads where the proposed order don't match fw and bw.
 *
 * @param nid node id to asses the happiness
 * @param order proposed order
 * @param trg_nt thread graph contains the threads generated using the reads
 * @return 3 long array with read count for happy, unhappy and disconnected reads [happy count, unhappy count, disconnected count]
*/
std::array<uint64_t,3> assess_node_happiness(sgNodeID_t nid, const std::unordered_map<sgNodeID_t , uint32_t> &order, const DistanceGraph& trg_nt);

/** @brief uses relative position propagation to create a total order for a connected component
 *
 * @param dg Graph containing the subcomponent to sort
 * @param cc set of ids of the component to sort ( like dg.get_connected_component() )
 * @return map of nodes to starting positions <node, position>
 * */
std::unordered_map<sgNodeID_t , int64_t > sort_cc(const DistanceGraph& dg, std::unordered_set<sgNodeID_t> cc);

/** @brief pop a node from a read thread
 *
 * dg is a threaded graph, the read thread is identified using the support id
 *
 * This function as is is not very fast, to get fast result use the queue and pop approach.
 *
 * @param dg graph to pop the node from
 * @param node_id node to pop
 * @param read thread to pop the node from
 * @return true if popped, false otherwise
 * */
bool pop_node(DistanceGraph& dg, sgNodeID_t node_id, uint64_t read);

/** @brief pop a node from all threads in the dg graph
 *  dg is a threaded graph, the read thread is identified using the support id
 *
 * @param dg graph to pop the node from
 * @param node_id node to pop
 * @return number of nodes popped
 */
int pop_node_from_all(DistanceGraph &dg, sgNodeID_t node_id);



/** @brief Takes a thread and pops all unhappy/disconnected nodes from it, returns an empty thread if theres too many of them
 *
 * @param thread
 * @param trg
 * @param max_unhappy
 * @param disconnection_rate
 * @return
 */
std::vector<NodePosition> make_thread_happy(const std::vector<NodePosition> &thread,const DistanceGraph & trg, int max_unhappy=1, float disconnection_rate=.3);
void make_all_threads_happy(LongReadsRecruiter & lrr, DistanceGraph &trg, int max_unhappy=1, float disconnection_rate=.3);

//class ThreadedGraphSorter {
//public:
//    // TODO read ids (rids) are uint64_t type always, should we define a rid type?? YES
//    ThreadedGraphSorter(const DistanceGraph & _dg){};
//
//    /*
//     * Returns a list of reads linking to or from the spec nodeview
//     * TODO: Change nv for node maybe!?
//     * */
////    std::vector<uint64_t > rids_from_node(NodeView nv);
////    uint64_t shared_reads(NodeView nv1, NodeView nv2);
//
////    bool pop_node(sgNodeID_t node_id, uint64_t read);
////    bool multipop_node(sgNodeID_t node_id, uint64_t read);
//
//    // Not implemented bc not used in the pipeline
////    DistanceGraph make_summary_connection_graph(DistanceGraph& mdg, int links_to_connect=3){};
//
//
//
//    const DistanceGraph & dg;
//};

enum HappyPosition{Nowhere,FrontFW,FrontBW,MiddleFW,MiddleBW,BackFW,BackBW};
class HappyInsertionSorter;//FW declaration
/**
 * this keeps a record of a Node's pre/post conditions and computes if the node is happy to be added to the order
 */
class NodeAdjacencies {
public:
    NodeAdjacencies (){};
    void mark_used(const sgNodeID_t & nid);
    void mark_unused(const sgNodeID_t & nid);

    /** @brief Evaluates if there is a happy position for the node
     *
     * The size of the prevs and nexts are compared to fw_prevs, bw_prevs,fw_nexts and bw_nexts=0. the counters are
     * incremented every time the central node is seen in the adjacency of another node, so a high concordance means
     * that most of the adjacency is shared between the node and the order.
     *
     * @param used_perc perc of the pre/next context that needs to be shared to consider a node happy
     * @return
     */
    HappyPosition happy_to_add(float used_perc); //TODO::expand to return a HappyPosition

    /** @brief Finds the last of the prevs and first of the nexts, signed for orientation
     *
     * @param sorter his
     * @return Position of the last prev and the first next
     */
    std::pair<int64_t,int64_t> find_happy_place(const HappyInsertionSorter& sorter);
    std::unordered_set<sgNodeID_t> prevs,nexts;
    uint64_t fw_prevs=0,bw_prevs=0,fw_nexts=0,bw_nexts=0;
};

/**
 * insertion sorter: needs to be founded by a relatively good set of ordered nodes (i.e. the local order of a node)
 * keeps a record of previous/next nodes to every node: if most previous/next nodes to a node are in the order, add the node to it
 */
class HappyInsertionSorter {
public:
    HappyInsertionSorter(ReadThreadsGraph& _rtg):rtg(_rtg){};
    /**@brief Computes each node adjacencies in the threads in an all v all manner (within the thread)
     * Fills adjacencies map, key is node value is a NodeAdjacency object with adj information
     *
     * @param min_links min number of links for a node to be considered adjacent
     * @param radius Max search radius to compute adjs (node count)
     */
    void compute_adjacencies(int min_links=2, int radius=20);

    /**
     * Adds node to node_positions as in position 0 (is updated later).
     * Updates the candidate pool
     * Marks the nodes as used in the adj list
     *
     * @param nid node to add
     */
    void add_node(sgNodeID_t nid); //figures out orientation, then adds with sign, then runs mark_used on all this node's NodeAdjacencies (only those will be affected), removes this node from candidates and adds its adjacencies

    /** @brief gets the adjacency object that describes the adj of the node
     *
     * @param nid central node
     * @return NodeAdjacencies object describing the adj of the node
     */
    NodeAdjacencies & get_node_adjacencies(sgNodeID_t nid);
    std::vector<sgNodeID_t> nodes_to_add(float used_perc=.9); //finds nodes to add, sorts them by inter-dependency (i.e. nodes coming first are not made happier by subsequent nodes).

    /** @brief Uses the node's happy place to put it in the order, shifts other nodes around if needed. Returns false if the node is unhappy.
     *
     * @param nid
     * @param used_perc
     * @param solve_floating_by_rtg
     * @return
     */
    bool insert_node(sgNodeID_t nid,float used_perc=0.9, int64_t at_position=0);
    //TODO: validate_order -> checks all nodes are still happy in their place.
    void reset_positions();
    void start_order_from_node(sgNodeID_t nid,float used_perc=0.9,bool cleanup_initial_order=true);
    bool start_order_from_list(std::vector<sgNodeID_t> nodes);
    void grow_order(float used_perc=0.9,uint64_t steps=UINT64_MAX);
    void grow_order2(float used_perc=0.9,uint64_t steps=UINT64_MAX, bool write_detailed_log=false);
    LocalOrder get_order() const;

    int64_t get_node_position(sgNodeID_t nid) const; //negative means reverse, but order is abs!!!!
    void remove_node_from_everywhere(sgNodeID_t nid);

    void dump_adjacencies(std::string filename);
    void load_adjacencies(std::string filename);


    LocalOrder local_order_from_node(sgNodeID_t nid,float perc=.9,bool cleanup_initial_order=true);

    ReadThreadsGraph& rtg;
    std::set<sgNodeID_t> candidates;
    std::unordered_map<sgNodeID_t,int64_t> node_positions;
    std::unordered_map<sgNodeID_t,NodeAdjacencies> adjacencies;
};

class TheGreedySorter {
public:
    TheGreedySorter(const DistanceGraph& _trg_nt, sgNodeID_t founding_node=0);
    void update_read_nodes_in_order();
    std::vector<uint64_t> node_belonging_scores(int64_t nid);
    std::vector<uint64_t > rids_from_node(NodeView nv);
    uint64_t shared_reads(NodeView nv1, NodeView nv2);
    bool pop_node(sgNodeID_t node_id, uint64_t read);
    bool remove_node(sgNodeID_t node_id);
    void start_from_read(uint64_t rid, int min_confirmation);
    std::pair<sgNodeID_t,sgNodeID_t> get_thread_ends(int64_t rid);
    std::pair<int, int> evaluate_read(uint64_t rid, bool print_pos);
    std::vector<int64_t> thread_nodes(int64_t rid);
    std::vector<int64_t> thread_node_positions(int64_t rid);
    std::vector<int32_t > evaluate_read_nodeorder(uint64_t rid, bool print_pos);
    void add_read(uint64_t rid, int min_confirmation=2);
//    void TheGreedySorter::extend_solution(int min_support=2, int min_shared=10, int min_new=10);

    void write_connected_nodes_graph(std::string filename);

    std::unordered_map<sgNodeID_t , int64_t > sort_graph();

    const DistanceGraph& trg_nt;
    DistanceGraph dg;

    std::vector<sgNodeID_t > all_nodes;
    std::unordered_set<uint64_t > all_reads;
    std::set<sgNodeID_t > used_nodes;
    std::vector<sgNodeID_t > used_reads;
    std::unordered_map<uint64_t, std::vector<sgNodeID_t >> read_ends;

private:
    bool order_is_valid=false;
    std::unordered_map<sgNodeID_t , int64_t > order;
    std::unordered_map<uint64_t,uint64_t> read_nodes_in_order;
};


/**
 * This class creates local orders for some of the nodes. Not every node will have its order generated
 */
class LocalOrderMaker{
public:
    LocalOrderMaker(ReadThreadsGraph &_rtg);

    void add_order(sgNodeID_t nid, const LocalOrder & o);
    int orders_in_region(const NodeAdjacencies &na);//median of # or orders in prevs/nexts
    void make_orders(int coverage, int min_links, int radius, float perc=.7, int min_adj = 10, int min_order_size=100);//multi-thread, 1 sorter per thread, creates orders until all nodes have coverage median in prevs/nexts
    std::map<sgNodeID_t,LocalOrder> local_orders;
    std::map<sgNodeID_t,std::vector<uint64_t>> node_local_orders;
    std::vector<uint64_t> node_local_orders_count;
    ReadThreadsGraph & rtg;
};

#endif //SDG_THREADEDGRAPHSORTER_H
