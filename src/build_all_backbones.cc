#include <iostream>
#include <sglib/workspace/WorkSpace.hpp>
#include <sglib/processors/LinkageUntangler.hpp>
#include <sglib/processors/HaplotypeConsensus.hpp>
#include "cxxopts.hpp"


int main(int argc, char **argv) {
    WorkSpace ws;

    std::string workspace_file;
    std::string output_prefix;
    uint32_t from(0), to(0);

    cxxopts::Options options("build_backbones", "Create consensus for backbones");
    try
    {
        options.add_options()
            ("help", "Print help")
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
        sglib::OutputLog() << "Building backbones from " << from << " to " << to << std::endl;

    ws.load_from_disk("nano10x_rg_lmpLR10x_mapped.bsgws");

    ws.long_read_mappers[0].read_filtered_mappings("fm_10K3.bsgfrm");
    ws.linked_read_mappers[0].read_tag_neighbours("tag_neighbours_1000_0.03.data");

    ws.long_read_mappers[0].update_indexes();
    ws.long_read_mappers[0].improve_filtered_mappings();

    sglib::OutputLog() << "Selecting backbones" << std::endl;

    auto u=LinkageUntangler(ws);
    auto mldg=u.make_longRead_multilinkage(ws.long_read_mappers[0]);
    u.select_multi_linkage_linear_anchors(mldg,5);
    auto ldg=u.make_nextselected_linkage(mldg);

    auto backbones=ldg.get_all_lines(2,10000);

    sglib::OutputLog() << "Done selecting backbones" << std::endl;

    ws.long_read_mappers[0].read_paths.resize(ws.long_read_mappers[0].filtered_read_mappings.size());

    if (to == 0) {
        to = backbones.size();
    }
    sglib::OutputLog() << "Building backbones from " << from << " to " << to << std::endl;
    for (uint32_t backbone=from; backbone <= to and backbone < backbones.size(); ++backbone) {
        sglib::OutputLog() << "Starting consensus for backbone " << backbone << std::endl;

        for (const auto &n: backbones[backbone]) {
            if (n != 0){
                std::cout << "seq"<<std::abs(n)<<",";
            }
        }
        std::cout << std::endl;

        HaplotypeConsensus haplotypeConsensus(ws, mldg, ldg, backbones[backbone]);

        bool read_paths_in_file=false;
        {
            std::ifstream orp_file("oriented_read_paths_backbone_"+std::to_string(backbone)+".orp");
            if (orp_file.good()){
                read_paths_in_file=true;
            }
        }

        if (!read_paths_in_file) {
            auto useful_read = ws.long_read_mappers[0].create_read_paths(backbones[backbone]);


            sglib::OutputLog() << "Created read paths for " << useful_read.size() << " reads" << std::endl;
            std::cout << "Read ids: ";
            std::copy(useful_read.cbegin(), useful_read.cend(), std::ostream_iterator<uint64_t>(std::cout, ", "));

            auto max_rid = std::max_element(useful_read.cbegin(), useful_read.cend());
            haplotypeConsensus.oriented_read_paths.resize(*max_rid + 1000);
            haplotypeConsensus.orient_read_paths(useful_read);
            sglib::OutputLog() << "Done orienting " << useful_read.size() << " read paths" << std::endl;

            haplotypeConsensus.write_read_paths("oriented_read_paths_backbone_" + std::to_string(backbone) + ".orp");

            // Clear temp structures
            ws.long_read_mappers[0].all_paths_between.clear();
            for (const auto &read: useful_read) {
                ws.long_read_mappers[0].read_paths[read].clear();
            }
            sglib::OutputLog() << "Done clearing structures" << std::endl;

        } else {
            haplotypeConsensus.read_read_paths("oriented_read_paths_backbone_"+std::to_string(backbone)+".orp");
        }

        haplotypeConsensus.build_line_path();

        sglib::OutputLog() << "Backbone " << backbone << " consensus: " << std::endl;
        for (const auto &n: haplotypeConsensus.backbone_filled_path) {
            if (n != 0){
                std::cout << "seq"<<std::abs(n)<<",";
            }
        }
        std::cout << std::endl;
        sglib::OutputLog() << "Done building line path" << std::endl;


        std::string consensus = haplotypeConsensus.consensus_sequence();
        std::ofstream backbone_consensus_fasta("consensus"+std::to_string(backbone)+".fasta");
        backbone_consensus_fasta << ">backbone_consensus_" << backbone << std::endl;
        backbone_consensus_fasta << consensus << std::endl;
        sglib::OutputLog() << "Done consensus for backbone " << backbone << std::endl;
    }

    return 0;
}