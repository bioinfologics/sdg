#include <iostream>
#include <fstream>
#include <sglib/WorkSpace.hpp>
#include <sglib/processors/Untangler.hpp>
#include <sglib/processors/FlowFollower.hpp>
#include "sglib/logger/OutputLog.h"
#include "cxxopts.hpp"

struct Counter
{
    struct value_type { template<typename T> value_type(const T&) { } };
    void push_back(const value_type&) { ++count; }
    size_t count = 0;
};

template<typename T1, typename T2>
size_t intersection_size(const T1& s1, const T2& s2)
{
    Counter c;
    set_intersection(s1.begin(), s1.end(), s2.begin(), s2.end(), std::back_inserter(c));
    return c.count;
}

int main(int argc, char * argv[]) {
    std::cout << "Welcome to bsg-untangler"<<std::endl<<std::endl;
    std::cout << "Git origin: " << GIT_ORIGIN_URL << " -> "  << GIT_BRANCH << std::endl;
    std::cout << "Git commit: " << GIT_COMMIT_HASH << std::endl<<std::endl;
    std::cout << "Executed command:"<<std::endl;
    for (auto i=0;i<argc;i++) std::cout<<argv[i]<<" ";
    std::cout<<std::endl<<std::endl;

    std::string workspace_file,output_prefix;
    sglib::OutputLogLevel=sglib::LogLevels::DEBUG;
    bool devel_code=false;
    try
    {
        cxxopts::Options options("bsg-untangler", "graph-based repeat resolution and haplotype separation");

        options.add_options()
                ("help", "Print help")
                ("w,workspace", "input workspace", cxxopts::value<std::string>(workspace_file))
                ("o,output", "output file prefix", cxxopts::value<std::string>(output_prefix))
                ("d,devel", "execute development code (different main!)", cxxopts::value<bool>(devel_code))
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

    std::cout<<std::endl;
    WorkSpace ws;
    sglib::OutputLog()<<"Loading Workspace..."<<std::endl;
    ws.load_from_disk(workspace_file);
    ws.add_log_entry("bsg-untangler run started");
    sglib::OutputLog()<<"Loading Workspace DONE"<<std::endl;
    if (ws.path_datastores.size()==0) {
        sglib::OutputLog()<<"Finishing early because there's no path_datastores[0]"<<std::endl;
        return 1;
    }

    if (!devel_code) {

        Untangler u(ws);

        u.analise_paths_through_nodes();
        FlowFollower ff(ws);

        ff.load_flows_from_ws();
        std::ofstream skf("skated_paths.fasta");
        for (auto sp:ff.skate_from_all(4, 10000)) {
            skf << ">" << sp.nodes[0] << "_" << sp.nodes.back() << std::endl << sp.get_sequence() << std::endl;
        }

        std::cout << "TODO: pop error-bp bubbles by kci and paths" << std::endl;
        u.pop_errors_by_ci_and_paths();
        std::cout
                << "TODO: Solve bubbly paths by kci (both local and bubly-path-total) and paths through collapses (expand without joins, update paths!)"
                << std::endl;
        std::cout << "TODO: Join unitigs and remap all reads" << std::endl;
        std::cout << "TODO: skate through crap, generate a copy of the solution and re-connect" << std::endl;
        std::cout << "TODO: Join unitigs and remap all reads" << std::endl;
        sglib::OutputLog() << "All DONE!!!" << std::endl;

        //u.analise_paths_through_nodes();
    }
    else {
        //==================== Development code (i.e. random tests!) ==============
        sglib::OutputLog()<<"Selecting nodes..."<<std::endl;
        std::vector<std::vector<std::pair<sgNodeID_t,uint32_t>>> neighbours;
        std::vector<std::set<bsg10xTag>> node_tags;
        neighbours.resize(ws.sg.nodes.size());
        node_tags.resize(ws.sg.nodes.size());
        uint64_t total_bp=0;
        auto nodes=ws.select_from_all_nodes(1000,1000000,20,200000, 0.5, 1.5);
        sglib::OutputLog()<<"Populating node tags..."<<std::endl;
        for (auto n:nodes) {
            total_bp+=ws.sg.nodes[n].sequence.size();
            for (auto t:ws.linked_read_mappers[0].get_node_tags(n)) node_tags[n].insert(t);
        }
        sglib::OutputLog()<<nodes.size()<<" selected totalling "<<total_bp<<"bp "<<std::endl;
        sglib::OutputLog()<<"Computing shared tags"<<std::endl;
#pragma omp parallel for shared(neighbours) schedule(static,50)
        for (auto i1=0; i1<nodes.size(); ++i1){
            auto n1=nodes[i1];
            for (auto i2=i1+1;i2<nodes.size();i2++){
                auto n2=nodes[i2];
                uint32_t shared=intersection_size(node_tags[n1],node_tags[n2]);
                if (shared>=4) {
#pragma omp critical
                    {
                    neighbours[n1].emplace_back(n2,shared);
                    neighbours[n2].emplace_back(n1,shared);
                    }
                }
            }
        }
        sglib::OutputLog()<<"Sorting shared tags"<<std::endl;
        for (auto &nn:neighbours){
            std::sort(nn.begin(),nn.end(),[]( const std::pair<sgNodeID_t,uint32_t> &a,
                                              const std::pair<sgNodeID_t,uint32_t> &b) { return a.second>b.second; });
        }
        sglib::OutputLog()<<"Dumping shared tags"<<std::endl;
        std::ofstream sto("shared_tags.txt");
        uint64_t linked=0;
        for (auto n:nodes){
            sto<<n<<" ("<<node_tags[n].size()<<")";
            if (neighbours[n].size()>0) ++linked;
            for (auto nn:neighbours[n]) sto<<", ["<<nn.first<<", "<<nn.second<<"]";
            sto<<std::endl;
            sto<<n<<">>>> seq"<<n;
            for (auto nn:neighbours[n]) sto<<", seq"<<nn.first;
            sto<<std::endl;


        }
        sglib::OutputLog()<<linked<<" nodes with neighbours"<<std::endl;

    }
    return 0;
}

