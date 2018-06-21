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
                    auto node = sgNodeID_t (graph_index->seq[regs0[j].rid].name);
                    LongReadMapping mapping = createMapping(readID, regs0, j, node);
                    thread_mappings[omp_get_thread_num()].emplace_back(mapping);
                }
            }
            if (n_regs0 > 1) {
                for (int j = 0; j < n_regs0 - 1; ++j) {
                    auto fromNode = sgNodeID_t (graph_index->seq[regs0[j].rid].name);
                    auto toNode = sgNodeID_t (graph_index->seq[regs0[j+1].rid].name);
                    LongReadMapping mapping = createMapping(readID, regs0, j, fromNode);
                    thread_mappings[omp_get_thread_num()].emplace_back(mapping);
                    mapping = createMapping(readID, regs0, j+1, toNode);
                    thread_mappings[omp_get_thread_num()].emplace_back(mapping);
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
    mapping.node = node * (regs0[j].rev == 0) ? 1 : -1;
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
