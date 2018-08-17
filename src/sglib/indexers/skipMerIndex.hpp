//
// Created by Luis Yanes (EI) on 16/03/2018.
//

#ifndef BSG_SKIPMERINDEX_HPP
#define BSG_SKIPMERINDEX_HPP


#include <sglib/graph/SequenceGraph.hpp>
#include <sglib/types/KmerTypes.hpp>
#include <sglib/factories/SkipMerIndexFactory.hpp>
#include <sglib/readers/FileReader.h>
#include <sglib/readers/SequenceGraphReader.h>
#include <cmath>

class SkipMerIndex {
    class kmerPosFactory : protected SkipMerFactory {
    public:
        kmerPosFactory(uint8_t k, uint8_t m, uint8_t n) : SkipMerFactory(k, m, n) {}

        void setFileRecord(FastaRecord &rec) {
            currentRecord = rec;
            for (auto ni=0; ni < N; ++ni) {
                fkmer[ni] =0;
                rkmer[ni] = 0;
                last_unknown[ni] = 0;
            }
            fi=0;
        }
        // TODO: Adjust for when K is larger than what fits in uint64_t!
        const bool next_element(std::vector<std::pair<uint64_t,graphPosition>> &mers){
            for (auto ni=0; ni < N; ++ni) {
                fkmer[ni] =0;
                rkmer[ni] = 0;
                last_unknown[ni] = 0;
            }
            fi=0;
            uint64_t p(0);
            graphPosition pos;
            while (p < currentRecord.seq.size()) {
                //fkmer: grows from the right (LSB)
                //rkmer: grows from the left (MSB)
                for (auto ni=0;ni<N;++ni) {
                    cycle_pos[ni]++;
                    if (cycle_pos[ni]==N)cycle_pos[ni]=0;
                    if (cycle_pos[ni]<M) {
                        fillKBuf(currentRecord.seq[p], fkmer[ni], rkmer[ni], last_unknown[ni]);
                    }
                    p++;
                    //if we are at p, the skip-mer started at p-S is now done

                    if (p < S-1) continue;

                    if (p == S - 1) fi = 0;
                    else {
                        ++fi;
                        if (fi == N) fi = 0;
                    }
                    if (unlikely(last_unknown[fi] + S > p)) continue;

                    if (fkmer[fi] <= rkmer[fi]) {
                        pos.node=currentRecord.id;
                        pos.pos=p;
                        mers.emplace_back(fkmer[fi],pos);
                    } else {
                        pos.node=-currentRecord.id;
                        pos.pos=p;
                        mers.emplace_back(rkmer[fi],pos);
                    }

                }
            }
            return false;
        }

    private:
        FastaRecord currentRecord;
        uint64_t bases;
    };

    using Map = std::unordered_map<uint64_t, graphStrandPos>;
    using const_iterator = std::unordered_map<uint64_t, graphStrandPos>::const_iterator;
    using pair = std::pair<uint64_t, graphStrandPos>;

    Map kmer_to_graphposition;
    uint8_t k=0;
    uint8_t m=0;
    uint8_t n=0;
    std::vector<uint64_t> unique_kmers_per_node;
    std::vector<uint64_t> total_kmers_per_node;

public:
    SkipMerIndex(){}
    SkipMerIndex(uint8_t k, uint8_t m, uint8_t n) :
            k(k), m(m), n(n)
    {
    }

    SkipMerIndex(const SequenceGraph& sg, uint8_t k, uint8_t m, uint8_t n) :
            k(k), m(m), n(n)
    {
        generate_index(sg);
    }

    void generate_index(const SequenceGraph & sg, bool verbose = true) {
        sglib::OutputLog() << "Generating index for " << sg.filename << std::endl;
        const std::string output_prefix("./");

        uint64_t total_k { 0 };
        total_kmers_per_node = std::vector<uint64_t>(sg.nodes.size(), 0);
        for (sgNodeID_t node = 0; node < sg.nodes.size(); node++) {
            auto sgnode = sg.nodes[node];
            if (sgnode.sequence.size() >= k) {
                auto n = sgnode.sequence.size() + 1 - k;
                total_k += n;
                total_kmers_per_node[node] = n;
            }
        }
        std::vector<std::pair<uint64_t,graphPosition>> kidxv;
        FastaRecord r;
        kidxv.reserve(total_k);
        kmerPosFactory kf( k, m, n );
        for (sgNodeID_t n = 1; n < sg.nodes.size(); ++n) {
            if (sg.nodes[n].sequence.size() >= k) {
                r.id = n;
                r.seq = sg.nodes[n].sequence;
                kf.setFileRecord(r);
                kf.next_element(kidxv);
            }
        }

        if (verbose) sglib::OutputLog(sglib::INFO)<<kidxv.size()<<" kmers in total"<<std::endl;
        if (verbose) sglib::OutputLog(sglib::INFO) << "  Sorting..."<<std::endl;
#ifdef _OPENMP
        __gnu_parallel::sort(kidxv.begin(),kidxv.end(),[](const std::pair<uint64_t,graphPosition> & a, const std::pair<uint64_t,graphPosition> & b){return a.first<b.first;});
#else
        std::sort(kidxv.begin(),kidxv.end(),[](const std::pair<uint64_t,graphPosition> & a, const std::pair<uint64_t,graphPosition> & b){return a.first<b.first;});
#endif

        sglib::OutputLog(sglib::INFO) << "  Merging..." << std::endl;
        if (verbose) sglib::OutputLog(sglib::INFO) << "  Merging..."<<std::endl;
        auto wi=kidxv.begin();
        auto ri=kidxv.begin();
        auto nri=kidxv.begin();
        while (ri<kidxv.end()){
            while (nri!=kidxv.end() and nri->first==ri->first) ++nri;
            if (nri-ri==1) {
                *wi=*ri;
                ++wi;
            }
            ri=nri;
        }

        kidxv.resize(wi - kidxv.begin());
        sglib::OutputLog(sglib::INFO) << kidxv.size() << " unique kmers in index, creating map" << std::endl;
        std::unordered_set<sgNodeID_t > seen_contigs;
        seen_contigs.reserve(sg.nodes.size());
        unique_kmers_per_node = std::vector<uint64_t>(sg.nodes.size(), 0);
        for (auto &kidx :kidxv) {
            kmer_to_graphposition[kidx.first] = { kidx.second.node, kidx.second.pos };
            unique_kmers_per_node[std::abs(kidx.second.node)] += 1;
            seen_contigs.insert(std::abs(kidx.second.node));
        }
        sglib::OutputLog(sglib::INFO) << seen_contigs.size() << " nodes with indexed kmers" <<std::endl;

        sglib::OutputLog() << "Done!" << std::endl;
    }

    const Map& getMap() {
        return kmer_to_graphposition;
    };

    const_iterator find(uint64_t hash) {
        return kmer_to_graphposition.find(hash);
    };

    const_iterator end() {
        return kmer_to_graphposition.cend();
    };

    void write_to_file(std::string &filename) {
        sglib::OutputLog() << "Dumping index" << std::endl;
        std::ofstream outf(filename, std::ios_base::binary);
        auto mapSize(kmer_to_graphposition.size());
        outf.write(reinterpret_cast<const char *>(&k), sizeof(k));
        outf.write(reinterpret_cast<const char *>(&m), sizeof(m));
        outf.write(reinterpret_cast<const char *>(&n), sizeof(n));
        outf.write(reinterpret_cast<const char *>(&mapSize), sizeof(mapSize));
        for (const auto &skm: kmer_to_graphposition){
            outf.write(reinterpret_cast<const char *>(&skm.first), sizeof(skm.first));
            outf.write(reinterpret_cast<const char *>(&skm.second), sizeof(skm.second));
        }
        sglib::OutputLog() << "Done!" << std::endl;
    }

    void load_from_file(std::string &filename) {
        sglib::OutputLog() << "Loading index" << std::endl;
        std::ifstream outf(filename, std::ios_base::binary);
        pair skm{0,graphStrandPos(0,0)};
        auto mapSize(kmer_to_graphposition.size());
        kmer_to_graphposition.reserve(mapSize);
        outf.read(reinterpret_cast<char *>(&k), sizeof(k));
        outf.read(reinterpret_cast<char *>(&m), sizeof(m));
        outf.read(reinterpret_cast<char *>(&n), sizeof(n));
        outf.read(reinterpret_cast<char *>(&mapSize), sizeof(mapSize));
        for (decltype(kmer_to_graphposition.size()) n = 0; n < mapSize; ++n){
            outf.read(reinterpret_cast<char *>(&skm.first), sizeof(skm.first));
            outf.read(reinterpret_cast<char *>(&skm.second), sizeof(skm.second));
            kmer_to_graphposition.emplace(skm.first, skm.second);
        }
        sglib::OutputLog() << "Done!" << std::endl;
    }


};


#endif //BSG_SKIPMERINDEX_HPP
