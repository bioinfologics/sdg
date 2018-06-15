//
// Created by Bernardo Clavijo (EI) on 11/02/2018.
//

#ifndef BSG_LINKEDREADMAPPER_HPP
#define BSG_LINKEDREADMAPPER_HPP

#include <map>

#include "sglib/graph/SequenceGraph.hpp"
#include <sglib/datastores/LinkedReadsDatastore.hpp>
#include <sglib/types/MappingTypes.hpp>

/**
 * @brief A mapper for linked reads from a LinkedReadsDatastore.
 *
 * Supports partial remapping of unmapped reads or of a selection list.
 */
class LinkedReadMapper {
    class StreamKmerFactory : public  KMerFactory {
    public:
        explicit StreamKmerFactory(uint8_t k) : KMerFactory(k){}
        inline void produce_all_kmers(const char * seq, std::vector<KmerIDX> &mers){
            // TODO: Adjust for when K is larger than what fits in uint64_t!
            last_unknown=0;
            fkmer=0;
            rkmer=0;
            auto s=seq;
            while (*s!='\0' and *s!='\n') {
                //fkmer: grows from the right (LSB)
                //rkmer: grows from the left (MSB)
                fillKBuf(*s, fkmer, rkmer, last_unknown);
                if (last_unknown >= K) {
                    if (fkmer <= rkmer) {
                        // Is fwd
                        mers.emplace_back(fkmer);
                        mers.back().contigID=1;
                    } else {
                        // Is bwd
                        mers.emplace_back(rkmer);
                        mers.back().contigID=-1;
                    }
                }
                ++s;
            }
        }
    };

public:
    LinkedReadMapper(SequenceGraph &_sg, LinkedReadsDatastore &_datastore) : sg(_sg),datastore(_datastore){
        reads_in_node.resize(sg.nodes.size());
    };
    void write(std::ofstream & output_file);
    void read(std::ifstream & input_file);
    void map_reads(std::unordered_set<uint64_t> const &  reads_to_remap={});
    void remap_all_reads();
    //void map_read(uint64_t readID);
    void remove_obsolete_mappings();
    /*void remap_reads();
    uint64_t process_reads_from_file(uint8_t, uint16_t, std::unordered_map<uint64_t , graphPosition> &, std::string , uint64_t, bool tags=false, std::unordered_set<uint64_t> const & reads_to_remap={});
    void save_to_disk(std::string filename);
    void load_from_disk(std::string filename);*/
    void print_stats(){};
    std::set<bsg10xTag> get_node_tags(sgNodeID_t n);

    SequenceGraph & sg;
    LinkedReadsDatastore &datastore;
    std::vector<std::vector<ReadMapping>> reads_in_node;
    std::vector<sgNodeID_t> read_to_node;//id of the main node if mapped, set to 0 to remap on next process

};


#endif //BSG_LINKEDREADMAPPER_HPP
