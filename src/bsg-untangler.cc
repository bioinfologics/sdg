#include <iostream>
#include <fstream>
#include <sglib/WorkSpace.hpp>
#include <sglib/GraphPartitioner.hpp>
#include <sglib/Scaffolder.hpp>
#include <sglib/processors/Untangler.hpp>
#include <sglib/processors/TagWalker.hpp>
#include "sglib/logger/OutputLog.h"
#include "cxxopts.hpp"


void walk_from(sgNodeID_t n, WorkSpace &ws){
    class StreamKmerFactory : public  KMerFactory {
    public:
        explicit StreamKmerFactory(uint8_t k) : KMerFactory(k){}
        inline void produce_all_kmers(const char * seq, std::unordered_set<uint64_t> &mers){
            // TODO: Adjust for when K is larger than what fits in uint64_t!
            last_unknown=0;
            fkmer=0;
            rkmer=0;
            auto s=seq;
            while (*s!='\0' and *s!='\n') {
                //fkmer: grows from the right (LSB)
                //rkmer: grows from the left (MSB)
                fillKBuf(*s, 0, fkmer, rkmer, last_unknown);
                if (last_unknown >= K) {
                    if (fkmer <= rkmer) {
                        // Is fwd
                        mers.insert(fkmer);
                    } else {
                        // Is bwd
                        mers.insert(rkmer);
                    }
                }
                ++s;
            }
        }
    };
    StreamKmerFactory skf(31);
    std::cout<<"Starting walk on "<<n<<"... "<<std::flush;
    auto tags=ws.linked_read_mappers[0].get_node_tags(n);
    std::cout<<tags.size()<<" tags... "<<std::flush;
    auto kmers=ws.linked_read_datastores[0].get_tags_kmers(31,3,tags);
    std::cout<<kmers.size()<<" kmers."<<std::endl;
    SequenceGraphPath p(ws.sg,{n});
    while (true){
        auto fwl=ws.sg.get_fw_links(p.nodes.back());
        if (fwl.empty()) break;
        sgNodeID_t best=0,second=0;
        double best_score=0,second_score=0;
        for (auto &l:fwl){
            std::unordered_set<uint64_t> lkmers, inters;
            skf.produce_all_kmers(ws.sg.nodes[(l.dest>0?l.dest:-l.dest)].sequence.c_str(),lkmers);
            std::set_intersection(lkmers.begin(),lkmers.end(),kmers.begin(),kmers.end(),std::inserter(inters,inters.end()));

            auto score=(double)(inters.size())/lkmers.size();
            std::cout<<"scoring transition to "<<l.dest<<": "<<inters.size()<<"/"<<lkmers.size()<<"="<<score<<std::endl;
            uint64_t hits=0;
            for (auto i:lkmers) if (kmers.count(i)) ++hits;
            score=(double)(hits)/lkmers.size();
            std::cout<<"scoring transition to "<<l.dest<<": "<<hits<<"/"<<lkmers.size()<<"="<<score<<std::endl;
            if (best==0 or score>best_score) {second=best;best=l.dest;second_score=best_score;best_score=score;}
            else if (second==0 or score>second_score) {second=l.dest;second_score=score;}
        }
        if (best_score==0) break;
        std::cout<<best_score<<" - "<<second_score<<std::endl;
        bool b=false;
        for (auto n:p.nodes) if (n==best) b=true;
        if (b) break;
        p.nodes.push_back(best);
        std::cout<<std::endl;
    }
    for (auto n:p.nodes) std::cout<<"seq"<<n<<", ";
}

int main(int argc, char * argv[]) {
    std::cout << "Welcome to bsg-untangler"<<std::endl<<std::endl;
    std::cout << "Git origin: " << GIT_ORIGIN_URL << " -> "  << GIT_BRANCH << std::endl;
    std::cout << "Git commit: " << GIT_COMMIT_HASH << std::endl<<std::endl;
    std::cout << "Executed command:"<<std::endl;
    for (auto i=0;i<argc;i++) std::cout<<argv[i]<<" ";
    std::cout<<std::endl<<std::endl;

    std::string workspace_file,output_prefix;
    bool stats_only=0,reindex_kci=0;
    uint64_t max_mem_gb=4;
    sglib::OutputLogLevel=sglib::LogLevels::DEBUG;
    bool repeats_first=1,pop_bubles=0,skip_remap=0,print_HSPNPs=1,haplotype_walk=0;
    std::string verbose_log="";
    try
    {
        cxxopts::Options options("bsg-untangler", "graph-based repeat resolution and haplotype separation");

        options.add_options()
                ("help", "Print help")
                ("w,workspace", "input workspace", cxxopts::value<std::string>(workspace_file))
                ("o,output", "output file prefix", cxxopts::value<std::string>(output_prefix));
                ("max_mem", "maximum_memory when mapping (GB, default: 4)", cxxopts::value<uint64_t>(max_mem_gb));

        options.add_options("Heuristics")
                ("skip_remap", "skip remapping reads on loading workspace (default false)", cxxopts::value<bool>(skip_remap))
                ("pop_short", "pop unsupported short bubbles first (default false)", cxxopts::value<bool>(pop_bubles))
                ("repeats_first", "solve_repeats before bubble phasing (default true)", cxxopts::value<bool>(repeats_first))
                ("print_hspnps", "print HSPNPs (default true)", cxxopts::value<bool>(print_HSPNPs))
                ("hwalk", "haplotype walk (default false)", cxxopts::value<bool>(haplotype_walk))
                ("heuristics_verbose_log", "dump heuristics verbose log to file", cxxopts::value<std::string>(verbose_log))
                ;



        auto result(options.parse(argc, argv));

        if (result.count("help"))
        {
            std::cout << options.help({"","Heuristics"}) << std::endl;
            exit(0);
        }

        if (result.count("w")!=1 or result.count("o")!=1) {
            throw cxxopts::OptionException(" please specify input workspace and output prefix");
        }



    } catch (const cxxopts::OptionException& e)
    {
        std::cout << "Error parsing options: " << e.what() << std::endl << std::endl
                <<"Use option --help to check command line arguments." << std::endl;
        exit(1);
    }


    WorkSpace ws;
    std::cout<<std::endl<<"=== Loading Workspace ==="<<std::endl;
    ws.load_from_disk(workspace_file);
    std::cout<<"Datastore size: "<<ws.linked_read_datastores[0].size()<<std::endl;
    ws.add_log_entry("bsg-untangler run started");
    ws.kci.compute_compression_stats();
    for (auto &m:ws.linked_read_mappers) m.memlimit=max_mem_gb*1024*1024*1024;
    if (skip_remap==false) {
        std::cout<<std::endl<<"=== Mapping reads ==="<<std::endl;
        for (auto &m:ws.linked_read_mappers) {
            m.update_graph_index();
            m.map_reads();
        }
        ws.add_log_entry("reads re-mapped to current graph");
        ws.dump_to_disk(output_prefix+"_mapped.bsgws");
    }




    std::cout<<std::endl<<"=== Scaffolding ==="<<std::endl;

    std::vector<PairedReadMapper> prmappers;
    std::ofstream verbose_log_file;
    if (verbose_log!=""){
        verbose_log_file.open(verbose_log);
    }
    /*if (pop_bubles) {
        scaff.pop_unsupported_shortbubbles();
        for (auto & lrm: scaff.lrmappers) {
            lrm.remove_obsolete_mappings();
            lrm.update_graph_index();
            lrm.map_reads();
        }
    }*/

    if (repeats_first) {
        std::cout << std::endl << "Solving trivial repeats" << std::endl;
        bool mod = true;
        int srpass = 0;
        while (mod) {
            Untangler u(ws);
            std::unordered_set<uint64_t> reads_to_remap;
            mod = false;
            ++srpass;
            uint64_t solved_paths=u.solve_canonical_repeats_by_tags(reads_to_remap);
            if (solved_paths>0) {
                mod=true;
                ws.kci.reindex_graph();
                //ws.sg.write_to_gfa(output_prefix + "_solved_repeats_" + std::to_string(srpass) + ".gfa",unsolved_repeats);
                std::cout << reads_to_remap.size() << " pairs of read to remap" << std::endl;
                ws.linked_read_mappers[0].update_graph_index();
                ws.linked_read_mappers[0].map_reads(reads_to_remap);
                ws.add_log_entry(std::to_string(solved_paths)+" new paths by solving canonical repeats, "+
                                         std::to_string(reads_to_remap.size()) +"reads re-mapped");
            }

        }
    }

    std::cout<<std::endl<<"Finding HSPNPs"<<std::endl;
    Untangler u(ws);
    auto hps=u.get_all_HSPNPs();
    if (print_HSPNPs) {
        std::cout << "Starting with " << hps.size() << " HSPNPs" << std::endl;
        std::ofstream hpsfile(output_prefix + "_hps.csv");
        auto hspnp_id = 1;
        for (auto hp:hps) { //TODO: print node numbers and alternative tags -> then cluster/find neighbours
            std::vector<prm10xTag_t> h1tags, h2tags;
            for (auto rm:ws.linked_read_mappers[0].reads_in_node[(hp.first > 0 ? hp.first : -hp.first)])
                h1tags.push_back(ws.linked_read_mappers[0].datastore.get_read_tag(rm.read_id));
            for (auto rm:ws.linked_read_mappers[0].reads_in_node[(hp.second > 0 ? hp.second : -hp.second)])
                h2tags.push_back(ws.linked_read_mappers[0].datastore.get_read_tag(rm.read_id));
            std::sort(h1tags.begin(), h1tags.end());
            std::sort(h2tags.begin(), h2tags.end());
            hpsfile << hspnp_id << ",A," << hp.first;
            for (auto t:h1tags) if (t) hpsfile << "," << t;
            hpsfile << std::endl;
            hpsfile << hspnp_id << ",B," << hp.second;
            for (auto t:h2tags) if (t) hpsfile << "," << t;
            hpsfile << std::endl;
            ++hspnp_id;
        }
    }
    if (haplotype_walk){
        for (auto hp:hps) {
            if (ws.sg.nodes[llabs(hp.first)].sequence.size()<2000 or ws.sg.nodes[llabs(hp.second)].sequence.size()<2000 ) continue;
            TagWalker tw(ws,hp);
            auto ct= tw.remove_crosstalk();
            if (ct>0) continue;
            tw.walk(.95,.1);
            //tw.dump_reads("HPSNP_"+std::to_string(llabs(hp.first))+"_"+std::to_string(llabs(hp.second)));


            //walk_from(hp.first,ws);
            //walk_from(hp.second,ws);
        }
    }




    ws.dump_to_disk(output_prefix+"_final.bsgws");
    ws.sg.write_to_gfa(output_prefix+"_scaffolded.gfa");
    return 0;
}

