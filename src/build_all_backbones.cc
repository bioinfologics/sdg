#include <iostream>
#include <sglib/workspace/WorkSpace.hpp>
#include <sglib/processors/LinkageUntangler.hpp>
#include <sglib/processors/HaplotypeConsensus.hpp>

int main(int argc, char **argv) {
    WorkSpace ws;

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


    sglib::OutputLog() << "Creating all paths between anchors" << std::endl;

//    for (const auto &b: backbones) {
//        uint32_t pos = 0;
//        auto &all_paths_between = ws.long_read_mappers[0].all_paths_between;
//        for (; pos < b.size()-1; pos++) {
//            const auto a1 = b[pos];
//            const auto a2 = b[pos+1];
//            auto tmp = ws.sg.find_all_paths_between(a1, a2, 20000, 40, false);
//            all_paths_between[std::make_pair(a1,a2)] = tmp;
//
//            for (auto &p:tmp){
//                p.reverse();
//            }
//
//            all_paths_between[std::make_pair(-a2,-a1)] = tmp;
//        }
//    }

    sglib::OutputLog() << "Done creating all paths between anchors " << std::endl;



    sglib::OutputLog() << "Begin pathing reads" << std::endl;

    ws.long_read_mappers[0].read_paths.resize(ws.long_read_mappers[0].filtered_read_mappings.size());
    ws.long_read_mappers[0].create_read_paths();

    sglib::OutputLog() << "Done pathing reads" << std::endl;

#pragma omp parallel for
    for (uint32_t backbone=0; backbone < backbones.size(); ++backbone) {
        sglib::OutputLog() << "Starting consensus for backbone " << backbone << std::endl;

        HaplotypeConsensus haplotypeConsensus(ws, mldg, ldg, backbones[backbone]);

        haplotypeConsensus.orient_read_paths();
        haplotypeConsensus.write_read_paths("oriented_read_paths_backbone_"+std::to_string(backbone)+".orp");
        haplotypeConsensus.build_line_path();
        std::string consensus = haplotypeConsensus.consensus_sequence();
        std::ofstream backbone_consensus_fasta("consensus"+std::to_string(backbone)+".fasta");
        backbone_consensus_fasta << ">backbone_consensus_" << backbone << std::endl;
        backbone_consensus_fasta << consensus << std::endl;
        sglib::OutputLog() << "Done consensus for backbone " << backbone << std::endl;
    }

    return 0;
}