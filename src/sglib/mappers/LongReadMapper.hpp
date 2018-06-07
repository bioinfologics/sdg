//
// Created by Luis Yanes (EI) on 12/02/2018.
//

#ifndef BSG_LONGREADMAPPER_H
#define BSG_LONGREADMAPPER_H


#include <iostream>
#include <sglib/factories/KMerIDXFactory.h>
#include <sglib/readers/FileReader.h>
#include <sglib/readers/SequenceGraphReader.h>
#include <sglib/SMR.h>
#include <sglib/factories/StrandedMinSketchFactory.h>
#include <sglib/indexers/minSketchIndex.hpp>
#include <sglib/datastores/LongReadsDatastore.hpp>
#include <sglib/types/MappingTypes.hpp>
#include <sglib/mappers/minimap2/minimap.h>


class LongReadMapper {
    SequenceGraph & sg;
    minSketchIndex index;
    LongReadsDatastore datastore;
    mm_idx_t *graph_index = nullptr;
    mm_mapopt_t opt;
    uint8_t k=15;
    uint8_t w=10;

    std::vector<LongReadMapping> mappings;
    std::vector< std::vector < std::vector<LongReadMapping>::size_type > > mappings_in_node;        /// Mappings matching node
    std::vector< std::vector < std::vector<LongReadMapping>::size_type > > read_to_mappings;    /// Nodes in the read, 0 or empty = unmapped

    void update_indexes_from_mappings();
    LongReadMapping createMapping(uint32_t readID, const mm_reg1_t *regs0, int j, long long int node) const;
    void printMatch(const mm_idx_t *mi, std::ofstream &matchOutput, uint32_t readID, const std::string &read_name,
                    int read_len, const mm_reg1_t *regs0, int j);

public:
    LongReadMapper(uint8_t k, uint8_t w, SequenceGraph &sg, LongReadsDatastore &ds) : sg(sg), k(k), w( ((w==0)? (uint8_t)(k*0.66f):w) ), index(sg, k, w), datastore(ds) {
        mappings_in_node.resize(sg.nodes.size());
        read_to_mappings.resize(datastore.size());

        std::vector<const char *> names(sg.nodes.size());
        std::vector<const char *> seqs(sg.nodes.size());
        for (std::vector<std::string>::size_type i = 1; i < seqs.size(); i++) {
            names[i] = sg.oldnames[i].data();
            seqs[i] = sg.nodes[i].sequence.data();
        }

        mm_mapopt_init(&opt);
        opt.flag |= MM_F_CIGAR;
    }

    void update_graph_index() {
        std::vector<const char *> names(sg.nodes.size());
        std::vector<const char *> seqs(sg.nodes.size());
        for (std::vector<std::string>::size_type i = 1; i < seqs.size(); i++) {
            names[i] = sg.oldnames[i].data();
            seqs[i] = sg.nodes[i].sequence.data();
        }
        if (graph_index != nullptr) {
            mm_idx_destroy(graph_index);
            graph_index = nullptr;
        }
        graph_index = mm_idx_str(w, k, 0, 14, static_cast<int>(seqs.size()-1), &seqs[1], &names[1]);
        mm_mapopt_update(&opt, graph_index);
    }

    ~LongReadMapper() {
        mm_idx_destroy(graph_index);
    }

    void map_reads(std::unordered_set<uint32_t> readIDs = {});

    void read_from_file(std::string filename) {
        // Read the mappings from file
        sglib::OutputLog() << "Reading long read mappings" << std::endl;
        std::ifstream outf(filename, std::ios_base::binary);
        auto mapSize(mappings.size());
        outf.read(reinterpret_cast<char *>(&k), sizeof(k));
        outf.read(reinterpret_cast<char *>(&w), sizeof(w));
        outf.read(reinterpret_cast<char *>(&mapSize), sizeof(mapSize));
        mappings.reserve(mapSize);
        outf.read(reinterpret_cast<char*>(mappings.data()), mappings.size()*sizeof(LongReadMapping));

        sglib::OutputLog() << "Updating read mapping indexes!" << std::endl;
        update_indexes_from_mappings();
        sglib::OutputLog() << "Done!" << std::endl;
    }
    void write_to_file(std::string filename) {
        // Write mappings to file
        sglib::OutputLog() << "Dumping long read mappings" << std::endl;
        std::ofstream outf(filename, std::ios_base::binary);
        auto mapSize(mappings.size());
        outf.write(reinterpret_cast<const char *>(&k), sizeof(k));
        outf.write(reinterpret_cast<const char *>(&w), sizeof(w));
        outf.write(reinterpret_cast<const char *>(&mapSize), sizeof(mapSize));
        outf.write(reinterpret_cast<const char*>(mappings.data()), mappings.size()*sizeof(LongReadMapping));
        sglib::OutputLog() << "Done!" << std::endl;
    }
};


#endif //BSG_LONGREADMAPPER_H
