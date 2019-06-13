#include <iostream>
#include <sdglib/workspace/WorkSpace.hpp>
#include <sdglib/processors/LinkageUntangler.hpp>
#include <sdglib/processors/HaplotypeConsensus.hpp>
#include "cxxopts.hpp"


int main(int argc, char **argv) {
    WorkSpace ws;

    std::string workspace_file;
    std::string mappings_file;
    std::string neighbours_file;
    bool use_checkpoints=false;
    std::string output_prefix;
    uint32_t from(0), to(0);

    cxxopts::Options options("build_backbones", "Create consensus for backbones");
    try
    {
        options.add_options()
            ("help", "Print help")
            ("w,workspace", "Load workspace", cxxopts::value(workspace_file))
            ("m,mappings", "Filtered mappings file", cxxopts::value(mappings_file))
            ("n,neighbours", "Pre-computed tag neighbours file", cxxopts::value(neighbours_file))
            ("c,checkpoints", "Enables the use of precomputed read paths if present", cxxopts::value(use_checkpoints)->default_value("false"))
            ("from", "backbone range start", cxxopts::value<uint32_t>(from))
            ("to", "backbone range end", cxxopts::value<uint32_t>(to));

        auto result(options.parse(argc, argv));

        if (result.count("help"))
        {
            std::cout << options.help({""}) << std::endl;
            exit(0);
        }
    } catch (const cxxopts::OptionException& e)
    {
        std::cout << e.what() << std::endl;
        std::cout << options.help({""}) << std::endl;
        exit(1);
    }

    if (from != 0 or to != 0)
        sdglib::OutputLog() << "Building backbones from " << from << " to " << to << std::endl;

    ws.load_from_disk(workspace_file);

    ws.long_read_mappers[0].read_filtered_mappings(mappings_file);
    ws.linked_read_mappers[0].read_tag_neighbours(neighbours_file);
    ws.long_read_mappers[0].update_indexes();

    sdglib::OutputLog() << "Selecting backbones" << std::endl;

    auto u=LinkageUntangler(ws);
    auto mldg=u.make_longRead_multilinkage(ws.long_read_mappers[0]);
    u.select_multi_linkage_linear_anchors(mldg,5);
    auto ldg=u.make_nextselected_linkage(mldg);

    auto backbones=ldg.get_all_lines(2,10000);

    sdglib::OutputLog() << "Done selecting backbones" << std::endl;

    if (to == 0) {
        to = backbones.size();
    }
    sdglib::OutputLog() << "Building backbones from " << from << " to " << to << std::endl;
    for (uint32_t backbone=from; backbone <= to and backbone < backbones.size(); ++backbone) {
        sdglib::OutputLog() << "Starting consensus for backbone " << backbone << std::endl;

        for (const auto &n: backbones[backbone]) {
            if (n != 0){
                std::cout << "seq"<<std::abs(n)<<",";
            }
        }
        std::cout << std::endl;

        bool read_paths_in_file=false;
        {
            std::ifstream orp_file("oriented_read_paths_backbone_"+std::to_string(backbone)+".orp");
            if (orp_file.good()){
                read_paths_in_file=true;
            }
        }

        ReadPathParams readPathParams;
        int min_votes=1;
        int min_path_nodes=2;
        int disconnected_distance=200;
        int min_distance=100;

        HaplotypeConsensus haplotypeConsensus(ws, mldg, ldg, backbones[backbone], readPathParams);
        std::string consensus;
        if (!read_paths_in_file or !use_checkpoints) {
            consensus = haplotypeConsensus.generate_consensus(min_votes, min_path_nodes, disconnected_distance, min_distance);

            haplotypeConsensus.write_read_paths("oriented_read_paths_backbone_" + std::to_string(backbone) + ".orp");

            // Clear temp structures
            ws.long_read_mappers[0].all_paths_between.clear();
            sdglib::OutputLog() << "Done clearing structures" << std::endl;

        } else {
            haplotypeConsensus.read_read_paths("oriented_read_paths_backbone_"+std::to_string(backbone)+".orp");
            haplotypeConsensus.build_line_path(min_votes, min_path_nodes);
            consensus = haplotypeConsensus.consensus_sequence(disconnected_distance, min_distance);
        }


        sdglib::OutputLog() << "Backbone " << backbone << " consensus: " << std::endl;
        for (const auto &n: haplotypeConsensus.backbone_filled_path) {
            if (n != 0){
                std::cout << "seq"<<std::abs(n)<<",";
            }
        }
        std::cout << std::endl;
        sdglib::OutputLog() << "Done building line path" << std::endl;


        std::ofstream backbone_consensus_fasta("consensus"+std::to_string(backbone)+".fasta");
        backbone_consensus_fasta << ">backbone_consensus_" << backbone << std::endl;
        backbone_consensus_fasta << consensus << std::endl;
        sdglib::OutputLog() << "Done consensus for backbone " << backbone << std::endl;
    }

    return 0;
}