//
// Created by Bernardo Clavijo (EI) on 24/02/2018.
//

#ifndef BSG_WORKSPACE_HPP
#define BSG_WORKSPACE_HPP


#include <sglib/datastores/LinkedReadsDatastore.hpp>
#include <sglib/mappers/LinkedReadMapper.hpp>
#include <sglib/datastores/PairedReadsDatastore.hpp>
#include <sglib/mappers/PairedReadMapper.hpp>
#include <sglib/datastores/LongReadsDatastore.hpp>
#include <sglib/mappers/LongReadMapper.hpp>
#include <sglib/datastores/PathsDatastore.hpp>
#include "sglib/graph/SequenceGraph.hpp"
#include "sglib/processors/KmerCompressionIndex.hpp"
#include <sglib/indexers/UniqueKmerIndex.hpp>

class LogEntry{
public:
    LogEntry(std::time_t t, std::string v, std::string tx):timestamp(t),bsg_version(std::move(v)),log_text(std::move(tx)){};
    std::time_t timestamp;
    std::string bsg_version;
    std::string log_text;
};

class WorkSpace {

public:
    WorkSpace() :
    kci(sg),
    uniqueKmerIndex(sg, 31),
    unique63merIndex(sg){};
    WorkSpace(const WorkSpace& that) = delete; //we definitely do not want copy constructors here, thank you

    /** @brief: Prints log content
     *
     */
    void print_log();

    /** @brief: Adds an entry to the log
     *
     * @param text text to be logged
     */
    void add_log_entry(std::string text);

    /** @brief: Dumps workspace to disk
     *
     * @param filename Name of the output file
     */
    void dump_to_disk(std::string filename);

    /** @brief: Loads workspace from disk
     *
     * @param filename Name of the workspace
     * @param log_only if true will only load the workspace log and return
     */
    void load_from_disk(std::string filename,bool log_only=false);

    //general operations
    /** @brief: Generates a unique kmer index at k=31 (see UniqueKmerIndex)
     * Is initialized with k31 in the workspace constructor
     * @param verbose
     */
    void create_index(bool verbose = true) { uniqueKmerIndex.generate_index(sg,verbose); }

    /** @brief: Generates a unique kmer index at k=63 (see Unique63merIndex)
     *
     * @param verbose
     */
    void create_63mer_index(bool verbose = true) { unique63merIndex.generate_index(sg,verbose); }

    /** @brief: Triggers the remap_all_reads for the pair read mappers and linked reads mapper of the workspace
     *  TODO: WARNING: This function will NOT remap the long reads
     */
    void remap_all();

    /** @brief: Triggers the remap all reads for the pair read mappers and linked reads mapper of the workspace usingn k=63
     * TODO: WARNING: This function will NOT remap the long reads
     */
    void remap_all63();
    //Projected operations with info from the graph

    /**@brief: Given a size range, a tag count range and a kci range returns a vector with the selected the nodeIDs from the nodes collection
     * if there is no linked_read_mappers the tags limit will be ignored
     * if there is no kci.graph_kmers collection the ci limits will be ignored
     *
     * TODO: for the linked reads only checks first mapped library
     * TODO: for the read counts only checks first counted kmer set
     * TODO: WARNING: this function does not check if a node is has been deleted, could potentially select non existng nodes.
     *
     * @param min_size Minimum node sequence length
     * @param max_size Maximum node sequence length
     * @param min_tags Minimum number of tags to consider a node
     * @param max_tags Maximum number of tags to consider a node
     * @param min_ci Minimum ci to consider a node
     * @param max_ci Maximum ci to consider a node
     * @return Vector of selected nodeIds
     */
    std::vector<sgNodeID_t>
    select_from_all_nodes(uint32_t min_size, uint32_t max_size, uint32_t min_tags, uint32_t max_tags, float min_ci, float max_ci);

    /** @brief: Returns a reference to the kci object (see KmerCompressionIndex)
     *
     * @return kci reference
     */
    KmerCompressionIndex& getKCI() {return kci;}

    /** @brief: Returns a reference to the graph object (see SequenceGraph)
     *
     * @return graph reference
     */
    SequenceGraph& getGraph() {return sg;}

    /** @brief: Vector storing LogEntry events
     *
     */
    std::vector<LogEntry> log;

    //All status classes are public, treat them with care anyway ;)
    /** @brief: graph object (see SequenceGraph)
     *
     */
    SequenceGraph sg;
    /** @brief: Unique Kmer index for the loaded graph (see UniqueKmerIndex)
     *
     */
    UniqueKmerIndex uniqueKmerIndex;
    /** @brief: Unique Kmer index at k=63 for the loaded graph (see Unique63merIndex)
     *
     */
    Unique63merIndex unique63merIndex;
    /** @brief: Paired read datastore collection (see PairedReadsDatastore)
     *
     */
    std::vector<PairedReadsDatastore> paired_read_datastores;
    /** @brief: Paired read mappers collection (see PairedReadMapper)
     *
     */
    std::vector<PairedReadMapper> paired_read_mappers;
    /** @brief: Linked read datastore collection (see LinkedReadsDatastore)
     *
     */
    std::vector<LinkedReadsDatastore> linked_read_datastores;
    /** @brief: Linked read mapper collection (see LinkedReadMapper)
     *
     */
    std::vector<LinkedReadMapper> linked_read_mappers;
    /** @brief: Long read datastore collection (see LongReadsDatastore)
     *
     */
    std::vector<LongReadsDatastore> long_read_datastores;
    /** @brief: Long read mapper collection (see LongReadMapper)
     *
     */
    std::vector<LongReadMapper> long_read_mappers;
    /** @brief: Collection of path datastores (see PathsDatastore)
     *
     */
    std::vector<PathsDatastore> path_datastores;
    /** @brief: Kmer compression index object (see KmerCompressionIndex)
     *
     */
    KmerCompressionIndex kci;
    std::string verbose_log="";

    static const bsgVersion_t min_compat;
    std::vector<LinkedReadMapper>& getLinkedReadMappers() {return linked_read_mappers;}
    std::vector<LinkedReadsDatastore>& getLinkedReadDatastores() {return linked_read_datastores;}
    std::vector<PairedReadMapper>& getPairedReadMappers() {return paired_read_mappers;}
    std::vector<PairedReadsDatastore>& getPairedReadDatastores() {return paired_read_datastores;}
    std::vector<LongReadMapper>& getLongReadMappers() {return long_read_mappers;}
    std::vector<LongReadsDatastore>& getLongReadDatastores() {return long_read_datastores;}
    std::vector<PathsDatastore>& getPathsDatastore() {return path_datastores;}
    std::vector<std::string> read_counts_header;
};


#endif //BSG_WORKSPACE_HPP
