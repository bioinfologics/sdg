//
// Created by Bernardo Clavijo (EI) on 24/02/2018.
//

#ifndef BSG_WORKSPACE_HPP
#define BSG_WORKSPACE_HPP


#include <sglib/datastores/LinkedReadsDatastore.hpp>
#include <sglib/mappers/LinkedReadMapper.hpp>
#include <sglib/datastores/PathsDatastore.hpp>
#include "SequenceGraph.h"
#include "KmerCompressionIndex.hpp"

class LogEntry{
public:
    LogEntry(std::time_t t, std::string v, std::string tx):timestamp(t),bsg_version(v),log_text(tx){};
    std::time_t timestamp;
    std::string bsg_version;
    std::string log_text;
};
class WorkSpace {

public:
    WorkSpace():kci(sg){};
    WorkSpace(const WorkSpace& that) = delete; //we definitely do not want copy constructors here, thank you
    void print_log();

    void add_log_entry(std::string text);

    void dump_to_disk(std::string filename);

    void load_from_disk(std::string filename,bool log_only=false);

    std::vector<LogEntry> log;

    //All status classes are public, treat them with care anyway ;)
    SequenceGraph sg;
    std::vector<LinkedReadsDatastore> linked_read_datastores;
    std::vector<LinkedReadMapper> linked_read_mappers;
    std::vector<PathsDatastore> path_datastores;
    KmerCompressionIndex kci;
    std::string verbose_log="";
};


#endif //BSG_WORKSPACE_HPP
