//
// Created by Bernardo Clavijo (EI) on 24/02/2018.
//

#pragma once

#include <sdglib/graph/SequenceDistanceGraph.hpp>

#include <sdglib/datastores/PairedReadsDatastore.hpp>
#include <sdglib/datastores/LinkedReadsDatastore.hpp>
#include <sdglib/datastores/LongReadsDatastore.hpp>

#include <sdglib/mappers/PairedReadsMapper.hpp>
#include <sdglib/mappers/LinkedReadsMapper.hpp>
#include <sdglib/mappers/LongReadsMapper.hpp>

#include <sdglib/indexers/UniqueKmerIndex.hpp>
#include <sdglib/workspace/Journal.hpp>
#include <sdglib/datastores/KmerCounter.hpp>

/**
 * The WorkSpace holds the information regarding a project in memory and can be written to a file for checkpointing purposes
 * It contains:
 *  - A SequenceDistanceGraph which holds all the sequences, and links between them. This is considered the base graph.
 *  - Multiple DistanceGraphs they contain alternative links between the sequences in the SDG nodes which can come from one or more sources of information (datastores)
 *  - Paired, Linked and Long reads datastores which contain all the sequencing information
 *  - KmerCounters which contain the kmers used for the count in the form of an index and multiple counts for those kmers. These counts can come from one or more data sources.
 *  - A Journal which contain details about any operation applied to the WorkSpace that transformed it in any way
 */
class WorkSpace {
public:
    explicit WorkSpace(const std::string & filename="");
    WorkSpace(const WorkSpace& that) = delete; //we definitely do not want copy constructors here, thank you
    void status();

    friend std::ostream& operator<<(std::ostream &os, const WorkSpace &ws);

    /**
     * @brief
     * Provides an overview of the information in the WorkSpace
     * @param level Base indentation level to use on the result
     * @param recursive Whether it should explore or not the rest of the hierarchy
     * @return
     * A text summary of the information contained in a WorkSpace
     */
    std::string ls(int level=0,bool recursive=true) const;

    /**
     * @brief Writes a disk version of the information held by that can be used in the future, it can be used as a form of checkpoint
     * @param filename Path to the file to write to disk
     */
    void dump_to_disk(std::string filename);

    /**
     * @brief
     * Load the workspace information from disk for further analysis
     * @param filename Path to the file to load from disk
     * @param log_only Only loads the Journal information
     */
    void load_from_disk(std::string filename,bool log_only=false);

    JournalOperation &add_operation(const std::string &name, const std::string &tool, const std::string &detail);


    /**
     * @brief
     * Add a new DistanceGraph and assign it a name
     * @param name Reference name to the DistanceGraph
     * @return Returns the DistanceGraph item added to the WorkSpace
     */
    DistanceGraph& add_distance_graph(const DistanceGraph &dg, const std::string &name);

    /**
     * @brief
     * Add a PairedReadsDatastore from an existing file and assign it a name
     * @param filename Path to the PairedReadsDatastore
     * @param name Reference name to the PairedReadsDatastore
     * @return Returns the PairedReadsDatastore item added to the WorkSpace
     */
    PairedReadsDatastore& add_paired_reads_datastore(const std::string &filename, const std::string &name="");

    /**
     * @brief
     * Add a LinkedReadsDatastore from an existing file and assign it a name
     * @param filename Path to the LinkedReadsDatastore
     * @param name Reference name to the LinkedReadsDatastore
     * @return Returns the LinkedReadsDatastore item added to the WorkSpace
     */
    LinkedReadsDatastore& add_linked_reads_datastore(const std::string &filename, const std::string &name="");

    /**
     * @brief
     * Add a LongReadsDatastore from an existing file and assign it a name
     * @param filename Path to the LongReadsDatastore
     * @param name Reference name to the LongReadsDatastore
     * @return Returns the LongReadsDatastore item added to the WorkSpace
     */
    LongReadsDatastore& add_long_reads_datastore(const std::string &filename, const std::string &name="");


    /**
     * @brief
     * Add a KmerCounter from an existing file and assign it a name
     * @param filename Path to the KmerCounter
     * @param name Reference name to the KmerCounter
     * @return Returns the KmerCounter item added to the WorkSpace
     */
    KmerCounter& add_kmer_counter(const std::string &filename,const std::string &name);

    /**
     * @brief
     * Add a new KmerCounter and assign it a name
     * @param name Reference name to the KmerCounter
     * @param k Size of the kmer used for the counts
     * @param count_mode Counts canonical or non-canonical kmers
     * @return Returns the KmerCounter item added to the WorkSpace
     */
    KmerCounter& add_kmer_counter(const std::string &name, const uint8_t k,
                                  KmerCountMode count_mode = KmerCountMode::Canonical);


    /**
     * @brief
     * Get a reference to a PairedReadsDatastore corresponding to the name
     * @param name Name of the PairedReadsDatastore
     * @return Returns the PairedReadsDatastore associated to the name
     */
    PairedReadsDatastore& get_paired_reads_datastore(const std::string &name);

    /**
     * @brief
     * Get a reference to a LinkedReadsDatastore corresponding to the name
     * @param name Name of the LinkedReadsDatastore
     * @return Returns the LinkedReadsDatastore associated to the name
     */
    LinkedReadsDatastore& get_linked_reads_datastore(const std::string &name);

    /**
     * @brief
     * Get a reference to a LongReadsDatastore corresponding to the name
     * @param name Name of the LongReadsDatastore
     * @return Returns the LongReadsDatastore associated to the name
     */
    LongReadsDatastore& get_long_reads_datastore(const std::string &name);

    /**
     * @brief
     * Get a reference to a DistanceGraph corresponding to the name
     * @param name Name of the DistanceGraph
     * @return Returns the DistanceGraph associated to the name
     */
    DistanceGraph& get_distance_graph(const std::string &name);

    /**
     * @brief
     * Get a reference to a KmerCounter corresponding to the name
     * @param name Name of the KmerCounter
     * @return Returns the KmerCounter associated to the name
     */
    KmerCounter& get_kmer_counter(const std::string &name);

    std::vector<std::string> list_distance_graphs();
    std::vector<std::string> list_paired_reads_datastores();
    std::vector<std::string> list_linked_reads_datastores();
    std::vector<std::string> list_long_reads_datastores();
    std::vector<std::string> list_kmer_counters();

    //Projected operations with info from the graph

    std::vector<sgNodeID_t>
    select_from_all_nodes(uint32_t min_size, uint32_t max_size, uint32_t min_tags, uint32_t max_tags, float min_ci, float max_ci);

    //All status classes are public, treat them with care anyway ;)
    SequenceDistanceGraph sdg;

    std::vector<PairedReadsDatastore> paired_reads_datastores;
    std::vector<LinkedReadsDatastore> linked_reads_datastores;
    std::vector<LongReadsDatastore> long_reads_datastores;
    std::vector<KmerCounter> kmer_counters;

    std::vector<DistanceGraph> distance_graphs;

    std::vector<JournalOperation> journal;

    static const sdgVersion_t min_compat;
};
