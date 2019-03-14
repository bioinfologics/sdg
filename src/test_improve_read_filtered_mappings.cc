#include <iostream>
#include <sglib/workspace/WorkSpace.hpp>
int main(int argc, char **argv) {

    WorkSpace ws;
    ws.load_from_disk("nano10x_rg_lmpLR10x_mapped.bsgws");

    ws.long_read_mappers[0].read_filtered_mappings("fm_10K3.bsgfrm");
    ws.long_read_mappers[0].update_indexes();
    ws.long_read_mappers[0].improve_read_filtered_mappings(100166);
    return 0;
}