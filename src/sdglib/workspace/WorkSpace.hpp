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
#include <sdglib/journal/OperationJournal.hpp>
#include <sdglib/datastores/KmerCounts.hpp>

class WorkSpace {

public:
    WorkSpace();
    explicit WorkSpace(const std::string & filename);
    WorkSpace(const WorkSpace& that) = delete; //we definitely do not want copy constructors here, thank you
    void status();


    void dump_to_disk(std::string filename);

    void dump_counts(std::string filename) {
        // Read until kmer counts, then append a new version of it
        std::fstream wsfile(filename);
        if (!wsfile.good()) {
            std::cerr << filename << " opening error: " << strerror(errno) << std::endl;
            throw std::runtime_error("Error opening " + filename);
        }

        //==== Check magic number ====
        sdgMagic_t magic;
        sdgVersion_t version;
        SDG_FILETYPE type;
        wsfile.seekg(sizeof(magic), std::ios_base::cur);
        wsfile.seekg(sizeof(version), std::ios_base::cur);
        wsfile.seekg(sizeof(type), std::ios_base::cur);

        uint64_t count,count2;
        //seek over operations
        wsfile.read((char *) &count, sizeof(count));
        for (auto i=0; i < count; i++) {
            sdglib::seek_string(wsfile); // Name
            sdglib::seek_string(wsfile); // Detail
            sdglib::seek_string(wsfile); // Tool
            wsfile.seekg(sizeof(std::time_t), std::ios_base::cur); // Timestamp
            wsfile.read((char *) &count2, sizeof(count2)); // Number of entries
            for (auto e = 0; e < count2; ++e){
                sdglib::seek_string(wsfile); // Entries
            }
        }

        // seek over graph
        wsfile.read((char *) &count,sizeof(count));
        for (auto i=0;i<count;++i) {
            wsfile.seekg(sizeof(NodeStatus)); // Status
            sdglib::seek_string(wsfile); // Sequence
        }
        Link link;
        sdglib::seek_flat_vectorvector(wsfile, link);

        // seek over distance graphs
        wsfile.read((char *) &count,sizeof(count));
        for (int i = 0; i < count; i++) {
            sdglib::seek_string(wsfile);
            sdglib::seek_flat_vectorvector(wsfile, link);
        }

        // seek over paired reads datastores
        wsfile.read((char *) &count,sizeof(count));
        for (int i = 0; i < count; i++) {
            sdglib::seek_string(wsfile);
            sdglib::seek_string(wsfile);

            wsfile.seekg(sizeof(magic), std::ios_base::cur);
            wsfile.seekg(sizeof(version), std::ios_base::cur);
            wsfile.seekg(sizeof(type), std::ios_base::cur);

            sgNodeID_t rtn;
            sdglib::seek_flat_vector(wsfile, rtn);
            ReadMapping rm;
            sdglib::seek_flat_vectorvector(wsfile, rm);
        }

        wsfile.read((char *) &count,sizeof(count));
        for (int i = 0; i < count; i++) {
            sdglib::seek_string(wsfile);
            sdglib::seek_string(wsfile);

            wsfile.seekg(sizeof(magic), std::ios_base::cur);
            wsfile.seekg(sizeof(version), std::ios_base::cur);
            wsfile.seekg(sizeof(type), std::ios_base::cur);

            sgNodeID_t rtn;
            sdglib::seek_flat_vector(wsfile, rtn);

            ReadMapping rin;
            sdglib::seek_flat_vectorvector(wsfile, rin);

            TagNeighbour tn;
            sdglib::seek_flat_vectorvector(wsfile, tn);
        }

        wsfile.read((char *) &count,sizeof(count));
        for (int i = 0; i < count; i++) {
            sdglib::seek_string(wsfile);
            sdglib::seek_string(wsfile);

            int8_t k;
            wsfile.seekg(sizeof(k), std::ios_base::cur);

            LongReadMapping mp;
            sdglib::seek_flat_vector(wsfile, mp);
        }

        count = kmer_counts.size();
        wsfile.write((char *) &count,sizeof(count));
        for (auto i = 0; i < count; i++) {
            kmer_counts[i].write(wsfile);
        }

    }

    void load_from_disk(std::string filename,bool log_only=false);

    OperationJournal &add_operation(const std::string &name, const std::string &tool, const std::string &detail);
    PairedReadsDatastore& add_paired_reads_datastore(const std::string &filename, const std::string &name="");
    LinkedReadsDatastore& add_linked_reads_datastore(const std::string &filename, const std::string &name="");
    LongReadsDatastore& add_long_reads_datastore(const std::string &filename, const std::string &name="");
    DistanceGraph& add_distance_graph(const DistanceGraph &dg, const std::string &name="");
    KmerCounts& add_kmer_counts_datastore(const std::string &name, const uint8_t k);


    PairedReadsDatastore& get_paired_reads_datastore(const std::string &name);
    LinkedReadsDatastore& get_linked_reads_datastore(const std::string &name);
    LongReadsDatastore& get_long_reads_datastore(const std::string &name);
    DistanceGraph& get_distance_graph(const std::string &name);
    KmerCounts& get_kmer_counts_datastore(const std::string &name);
    OperationJournal& get_operation(const std::string &name);

    //general operations

    void remap_all();
    void remap_all63();
    //Projected operations with info from the graph

    std::vector<sgNodeID_t>
    select_from_all_nodes(uint32_t min_size, uint32_t max_size, uint32_t min_tags, uint32_t max_tags, float min_ci, float max_ci);

    //All status classes are public, treat them with care anyway ;)
    SequenceDistanceGraph sdg;

    std::vector<PairedReadsDatastore> paired_reads_datastores;
    std::vector<LinkedReadsDatastore> linked_reads_datastores;
    std::vector<LongReadsDatastore> long_reads_datastores;
    std::vector<KmerCounts> kmer_counts;

    std::vector<DistanceGraph> distance_graphs;

    std::vector<OperationJournal> operation_journals;

    static const sdgVersion_t min_compat;
};
