//
// Created by Luis Yanes (EI) on 22/03/2018.
//

#include <sglib/mappers/LongReadMapper.hpp>
#include <sglib/utilities/omp_safe.hpp>

void LongReadMapper::map_reads(std::unordered_set<uint32_t> readIDs) {
    std::vector<std::vector<LongReadMapping>> thread_mappings(omp_get_max_threads());
#pragma omp parallel
    {
        mm_tbuf_t *buf = mm_tbuf_init();
#pragma omp for
        for (uint32_t readID = 1; readID < datastore.size(); ++readID) {
            if (!((readIDs.size()>0 and readIDs.count(readID)>0) or (readIDs.empty() and read_to_mappings[readID].empty())))
                continue;

            std::string read_seq(datastore.get_read_sequence(readID));
            std::string read_name(std::to_string(readID));
            int read_len(static_cast<int>(read_seq.length()));
            int n_regs0;
            mm_reg1_t *regs0 = mm_map(graph_index, read_len, read_seq.data(), &n_regs0, buf, &opt, read_name.data());

            if (n_regs0<=1) {
                for (int j = 0; j < n_regs0; j++) {
                    auto node = regs0[j].rid + 1;
                    LongReadMapping mapping = createMapping(readID, regs0, j, node);
                    thread_mappings[omp_get_thread_num()].emplace_back(mapping);
                }
            }
            if (n_regs0 > 1) {
                for (int j = 0; j < n_regs0 - 1; ++j) {
                    auto fromNode = regs0[j].rid+1;
                    auto toNode =   regs0[j+1].rid+1;
                    LongReadMapping fromMapping = createMapping(readID, regs0, j, fromNode);
                    LongReadMapping toMapping = createMapping(readID, regs0, j+1, toNode);

                    if (link_is_valid(fromMapping, toMapping)) {
                        thread_mappings[omp_get_thread_num()].emplace_back(fromMapping);
                        thread_mappings[omp_get_thread_num()].emplace_back(toMapping);
                    }
                }
            }
            for (int i = 0; i<n_regs0;i++) free(regs0[i].p);
            free(regs0);
        }
        mm_tbuf_destroy(buf);
    }
    for (int thread = 0; thread<omp_get_max_threads(); thread++) {
        mappings.reserve(thread_mappings.size());
        for (const auto &mapping : thread_mappings[thread]) {
            mappings.emplace_back(mapping);
        }
    }

    update_indexes_from_mappings();
}

LongReadMapping LongReadMapper::createMapping(uint32_t readID, const mm_reg1_t *regs0, int j, long long int node) const {
    LongReadMapping mapping;
    mapping.node = node * ((regs0[j].rev == 0) ? 1 : -1);
    mapping.nStart = regs0[j].rs;
    mapping.nEnd = regs0[j].re;
    mapping.qStart = regs0[j].qs;
    mapping.qEnd = regs0[j].qe;
    mapping.read_id = readID;
    mapping.score = regs0[j].score;
    return mapping;
}

void LongReadMapper::printMatch(const mm_idx_t *mi, std::ofstream &matchOutput, uint32_t readID,
                                const std::string &read_name, int read_len, const mm_reg1_t *regs0, int j) {
    matchOutput << readID
                << "\t" << read_name
                << "\t" << read_len
                << "\t" << j + 1
                << "\t" << regs0[j].rid
                << "\t" << mi->seq[regs0[j].rid].name
                << "\t" << mi->seq[regs0[j].rid].len
                << "\t" << regs0[j].rs
                << "\t" << regs0[j].re
                << "\t" << regs0[j].qs
                << "\t" << regs0[j].qe
                << "\t" << "+-"[regs0[j].rev]
                << "\t" << regs0[j].mlen
                << "\t" << regs0[j].blen
                << "\t" << regs0[j].mapq
                << "\n";
}

void LongReadMapper::update_indexes_from_mappings() {
    for (std::vector<LongReadMapping>::const_iterator mappingItr = mappings.cbegin(); mappingItr != mappings.cend(); ++mappingItr) {
        auto index = (unsigned long)std::distance(mappings.cbegin(), mappingItr);
        read_to_mappings[mappingItr->read_id].push_back(index);
        mappings_in_node[std::abs(mappingItr->node)].push_back(index);
    }
}

LongReadMapper::LongReadMapper(SequenceGraph &sg, LongReadsDatastore &ds, uint8_t k, uint8_t w)
        : sg(sg), k(k), w(((w == 0) ? (uint8_t)(k * 0.66f) : w) ), datastore(ds) {
    mappings_in_node.resize(sg.nodes.size());
    read_to_mappings.resize(datastore.size());

    mm_mapopt_init(&opt);
    opt.flag |= MM_F_CIGAR;
}

void LongReadMapper::update_graph_index() {
    std::vector<std::string> names(sg.nodes.size());
    std::vector<const char *> pt_names(sg.nodes.size());
    std::vector<const char *> seqs(sg.nodes.size());
    for (std::vector<std::string>::size_type i = 1; i < seqs.size(); i++) {
        names[i] = std::to_string(i);
        pt_names[i] = names[i].data();
        seqs[i] = sg.nodes[i].sequence.data();
    }
    if (graph_index != nullptr) {
        mm_idx_destroy(graph_index);
        graph_index = nullptr;
    }
    graph_index = mm_idx_str(w, k, 0, 14, static_cast<int>(seqs.size()-1), &seqs[1], &pt_names[1]);
    mm_mapopt_update(&opt, graph_index);
}

LongReadMapper::~LongReadMapper() {
    mm_idx_destroy(graph_index);
}

bool LongReadMapper::link_is_valid(const LongReadMapping &fromMapping, const LongReadMapping &toMapping) {
    // Check that the mappings are larger than 200bp
    if (fromMapping.qEnd-fromMapping.qStart < 80 or toMapping.qEnd-toMapping.qStart < 80) {
        return false;
    }
    // Check the query is not overlapped by more than 60%
    if (toMapping.qStart <= fromMapping.qEnd) {
        auto olpSize = std::max(fromMapping.qStart,toMapping.qStart) - std::min(fromMapping.qEnd,toMapping.qEnd);
        if (olpSize > std::abs(fromMapping.qEnd - fromMapping.qStart)*.8 or olpSize > std::abs(toMapping.qEnd-toMapping.qStart)*.8) {
            return false;
        }
    }

    // Check that the score is over 90% on both ends
//        if (fromMapping.score < std::abs(fromMapping.qEnd-fromMapping.qStart)*.9f or
//            toMapping.score < std::abs(toMapping.qEnd - toMapping.qStart)*.9f) {
//            return false;
//        }

    return true;
}

void LongReadMapper::read(std::string filename) {
    // Read the mappings from file
    sglib::OutputLog() << "Reading long read mappings" << std::endl;
    std::ifstream inf(filename, std::ios_base::binary);
    auto mapSize(mappings.size());
    inf.read(reinterpret_cast<char *>(&k), sizeof(k));
    inf.read(reinterpret_cast<char *>(&w), sizeof(w));
    inf.read(reinterpret_cast<char *>(&mapSize), sizeof(mapSize));
    mappings.reserve(mapSize);
    inf.read(reinterpret_cast<char*>(mappings.data()), mappings.size()*sizeof(LongReadMapping));

    sglib::OutputLog() << "Updating read mapping indexes!" << std::endl;
    update_indexes_from_mappings();
    sglib::OutputLog() << "Done!" << std::endl;
}

void LongReadMapper::read(std::ifstream &inf) {
    auto mapSize(mappings.size());
    inf.read(reinterpret_cast<char *>(&k), sizeof(k));
    inf.read(reinterpret_cast<char *>(&w), sizeof(w));
    inf.read(reinterpret_cast<char *>(&mapSize), sizeof(mapSize));
    mappings.resize(mapSize);
    inf.read(reinterpret_cast<char*>(mappings.data()), mappings.size()*sizeof(LongReadMapping));

    sglib::OutputLog() << "Updating read mapping indexes!" << std::endl;
    update_indexes_from_mappings();
    sglib::OutputLog() << "Done!" << std::endl;
}

void LongReadMapper::write(std::string filename) {
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

void LongReadMapper::write(std::ofstream &ofs) {
    auto mapSize(mappings.size());
    ofs.write(reinterpret_cast<char *>(&k), sizeof(k));
    ofs.write(reinterpret_cast<char *>(&w), sizeof(w));
    ofs.write(reinterpret_cast<char *>(&mapSize), sizeof(mapSize));
    mappings.reserve(mapSize);
    ofs.write(reinterpret_cast<char*>(mappings.data()), mappings.size()*sizeof(LongReadMapping));
    sglib::OutputLog() << "Done!" << std::endl;
}

LongReadMapper LongReadMapper::operator=(const LongReadMapper &other) {
    if (&sg != &other.sg and &datastore != &other.datastore) { throw ("Can only copy paths from the same SequenceGraph"); }
    if (&other == this) {
        return *this;
    }
    mappings = other.mappings;
    update_indexes_from_mappings();
    return *this;
}
