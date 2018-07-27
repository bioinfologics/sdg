//
// Created by Luis Yanes (EI) on 22/03/2018.
//

#include <sglib/mappers/LongReadMapper.hpp>
#include <sglib/utilities/omp_safe.hpp>

void LongReadMapper::map_reads(std::unordered_set<uint32_t> readIDs) {
    std::vector<std::vector<LongReadMapping>> thread_mappings(omp_get_max_threads());
#pragma omp parallel
    {
        ASMKMerFactory skf(17);
        std::vector<kmerPos > read_kmers;
#pragma omp for
        for (uint32_t readID = 1; readID < datastore.size(); ++readID) {
            if (!((readIDs.size()>0 and readIDs.count(readID)>0) or (readIDs.empty() and read_to_mappings[readID].empty())))
                continue;

            std::string read_seq(datastore.get_read_sequence(readID));
            std::string read_name(std::to_string(readID));
            int read_len(static_cast<int>(read_seq.length()));

            int window_size = 500;
            int window_slide = window_size/4;

            std::cout << "Processing Read ID " << readID << " " << read_seq.size() << std::endl;
            std::map<uint32_t , uint32_t > node_score;
            //===== Create a vector of nodes for every k-mer in the read.
            read_kmers.clear();
            skf.create_kmers(read_seq, 0, read_kmers); //XXX: Ns will produce "deletion-windows"
            std::cout<<"creating the vector of vectors of node matches... SLOOOOOOW"<<std::endl;
            std::vector<std::vector<uint32_t>> node_matches(read_kmers.size());
            for (auto i=0;i<read_kmers.size();++i){
                auto first = std::lower_bound(assembly_kmers.begin(), assembly_kmers.end(), read_kmers[i].kmer);
                for (auto it = first; it != assembly_kmers.end() && it->kmer == read_kmers[i].kmer; ++it) {
                    node_matches[i].emplace_back(it->contigID);
                    ++node_score[it->contigID];
                }
            }
            std::cout<<"DONE!"<<std::endl;



            //===== Rank nodes by how many kmer-hits they have on the read, choose the top N
            std::cout<<"Ranking nodes"<<std::endl;
            std::vector<std::pair<uint32_t , uint32_t>> score_node;
            for (auto &ns:node_score){
                score_node.emplace_back(ns.second,ns.first);
            }
            std::sort(score_node.rbegin(),score_node.rend());
            if (score_node.size()>100) score_node.resize(100);
            std::cout<<"DONE"<<std::endl;


            //===== For every window (sliding by win_size/2) save total of kmers in every node of the top N
            std::cout<<"Collecting votes for the winner nodes"<<std::endl;
            std::ofstream readout("read" + std::to_string(readID/2) + ".csv");
            readout<<"window#";
            for (auto n:score_node) readout<<",ctg"<<n.second;
            readout<<std::endl;

            //auto num_windows((read_kmers.size()-window_size)/overlap + 1);

            for (int w_i = 0; w_i * window_slide + window_size < read_kmers.size(); ++w_i) {
                std::map<uint32_t,uint32_t> win_node_score;
                auto win_start(w_i * window_slide);
                auto win_end(w_i * window_slide + window_size);

                for (int wp_i = win_start; wp_i < win_end and wp_i < read_kmers.size(); wp_i++) {
                    for (auto matchnode:node_matches[wp_i]) ++win_node_score[matchnode];
                }
                std::cout << "Window " << w_i+1 << " " << win_start << " to " << win_end << std::endl;
                readout<<w_i;
                for (auto n:score_node) readout<<","<<win_node_score[n.second];
                readout<<std::endl;
            }
            std::cout<<"DONE"<<std::endl;

//            if (n_regs0<=1) {
//                for (int j = 0; j < n_regs0; j++) {
//                    auto node = regs0[j].rid + 1;
//                    LongReadMapping mapping = createMapping(readID, regs0, j, node);
//                    thread_mappings[omp_get_thread_num()].emplace_back(mapping);
//                }
//            }
//            if (n_regs0 > 1) {
//                for (int j = 0; j < n_regs0 - 1; ++j) {
//                    auto fromNode = regs0[j].rid+1;
//                    auto toNode =   regs0[j+1].rid+1;
//                    LongReadMapping fromMapping = createMapping(readID, regs0, j, fromNode);
//                    LongReadMapping toMapping = createMapping(readID, regs0, j+1, toNode);
//
//                    if (link_is_valid(fromMapping, toMapping)) {
//                        thread_mappings[omp_get_thread_num()].emplace_back(fromMapping);
//                        thread_mappings[omp_get_thread_num()].emplace_back(toMapping);
//                    }
//                }
//            }
//            for (int i = 0; i<n_regs0;i++) free(regs0[i].p);
//            free(regs0);
        }
//        mm_tbuf_destroy(buf);

    }
    for (int thread = 0; thread<omp_get_max_threads(); thread++) {
        mappings.reserve(thread_mappings.size());
        for (const auto &mapping : thread_mappings[thread]) {
            mappings.emplace_back(mapping);
        }
    }

    update_indexes_from_mappings();
}

void LongReadMapper::update_indexes_from_mappings() {
    for (std::vector<LongReadMapping>::const_iterator mappingItr = mappings.cbegin(); mappingItr != mappings.cend(); ++mappingItr) {
        auto index = (unsigned long)std::distance(mappings.cbegin(), mappingItr);
        read_to_mappings[mappingItr->read_id].push_back(index);
        mappings_in_node[std::abs(mappingItr->node)].push_back(index);
    }
}

LongReadMapper::LongReadMapper(SequenceGraph &sg, LongReadsDatastore &ds, uint8_t k)
        : sg(sg), k(k), datastore(ds) {
    mappings_in_node.resize(sg.nodes.size());
    read_to_mappings.resize(datastore.size());
}

void LongReadMapper::update_graph_index() {
    assembly_kmers.reserve(110000000);

    ASMKMerFactory skf(17);
    std::vector<kmerPos > contig_kmers;

    for (sgNodeID_t n = 1; n < sg.nodes.size(); ++n) {
        if (sg.nodes[n].sequence.size() >= k) {
            contig_kmers.clear();
            skf.create_kmers(sg.nodes[n].sequence, n, contig_kmers);
            for (const auto &kmer:contig_kmers) {
                assembly_kmers.emplace_back(kmer.kmer, kmer.contigID, kmer.offset);
            }
        }
    }
    std::sort(assembly_kmers.begin(),assembly_kmers.end(), kmerPos::byKmerContigOffset());

}

LongReadMapper::~LongReadMapper() {}

void LongReadMapper::read(std::string filename) {
    // Read the mappings from file
    sglib::OutputLog() << "Reading long read mappings" << std::endl;
    std::ifstream inf(filename, std::ios_base::binary);
    auto mapSize(mappings.size());
    inf.read(reinterpret_cast<char *>(&k), sizeof(k));
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
    outf.write(reinterpret_cast<const char *>(&mapSize), sizeof(mapSize));
    outf.write(reinterpret_cast<const char*>(mappings.data()), mappings.size()*sizeof(LongReadMapping));
    sglib::OutputLog() << "Done!" << std::endl;
}

void LongReadMapper::write(std::ofstream &ofs) {
    auto mapSize(mappings.size());
    ofs.write(reinterpret_cast<char *>(&k), sizeof(k));
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
