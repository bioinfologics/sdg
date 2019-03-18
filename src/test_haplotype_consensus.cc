#include <iostream>
#include <sglib/workspace/WorkSpace.hpp>
#include <sglib/processors/LinkageUntangler.hpp>
#include <sglib/processors/HaplotypeConsensus.hpp>

int main(int argc, char **argv) {
    std::cout << "Hello world" << std::endl;

    WorkSpace ws;
    ws.load_from_disk("nano10x_rg_lmpLR10x_mapped.bsgws");

    ws.long_read_mappers[0].read_filtered_mappings("fm_10K3.bsgfrm");
    ws.long_read_mappers[0].update_indexes();
    ws.long_read_mappers[0].improve_filtered_mappings();

    ws.linked_read_mappers[0].read_tag_neighbours("tag_neighbours_1000_0.03.data");


    auto u = LinkageUntangler(ws);
    auto imldg =u.make_longRead_multilinkage(ws.long_read_mappers[0]);
    u.select_multi_linkage_linear_anchors(imldg,5);

    auto ildg1 = u.make_nextselected_linkage(imldg); // This should be the LDG of the backbones

    auto ilines = ildg1.get_all_lines(2,10000); // Confirm these are backbones?!


    for (int i = 0; i < 10; i++){
        HaplotypeConsensus haplotypeConsensus(ws, imldg, ildg1, ilines[i]);
        haplotypeConsensus.use_long_reads_from_file("reads_in_iline_"+std::to_string(i)+".fasta");
        haplotypeConsensus.orient_read_paths();
        haplotypeConsensus.build_line_path();
        haplotypeConsensus.consensus_sequence();
    }
    return 0;
}