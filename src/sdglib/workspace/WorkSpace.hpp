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

#include <sdglib/processors/KmerCompressionIndex.hpp>
#include <sdglib/indexers/UniqueKmerIndex.hpp>


class LogEntry{
public:
    LogEntry(std::time_t t, std::string v, std::string tx):timestamp(t),bsg_version(std::move(v)),log_text(std::move(tx)){};
    std::time_t timestamp;
    std::string bsg_version;
    std::string log_text;
};

class WorkSpace {

public:
    WorkSpace() : kci(sdg){}
    explicit WorkSpace(const std::string & filename) : kci(sdg){
        load_from_disk(filename);
    }
    WorkSpace(const WorkSpace& that) = delete; //we definitely do not want copy constructors here, thank you
    void print_log();

    void add_log_entry(std::string text);

    void dump_to_disk(std::string filename);

    void load_from_disk(std::string filename,bool log_only=false);

    PairedReadsDatastore& add_paired_reads_datastore(const std::string &name, const std::string &filename);
    LinkedReadsDatastore& add_linked_reads_datastore(const std::string &name, const std::string &filename);
    LongReadsDatastore& add_long_reads_datastore(const std::string &name, const std::string &filename);
    DistanceGraph& add_distance_graph(const std::string &name, const DistanceGraph &dg);


    PairedReadsDatastore& get_paired_reads_datastore(const std::string &name);
    LinkedReadsDatastore& get_linked_reads_datastore(const std::string &name);
    LongReadsDatastore& get_long_reads_datastore(const std::string &name);
    DistanceGraph& get_distance_graph(const std::string &name);

    //general operations

    void remap_all();
    void remap_all63();
    //Projected operations with info from the graph

    std::vector<sgNodeID_t>
    select_from_all_nodes(uint32_t min_size, uint32_t max_size, uint32_t min_tags, uint32_t max_tags, float min_ci, float max_ci);

    std::vector<LogEntry> log;

    //All status classes are public, treat them with care anyway ;)
    SequenceDistanceGraph sdg;

    std::vector<PairedReadsDatastore> paired_read_datastores;
    std::vector<LinkedReadsDatastore> linked_read_datastores;
    std::vector<LongReadsDatastore> long_read_datastores;

    std::vector<DistanceGraph> distance_graphs;

    KmerCompressionIndex kci;

    static const sdgVersion_t min_compat;
    std::vector<std::string> read_counts_header;
};
