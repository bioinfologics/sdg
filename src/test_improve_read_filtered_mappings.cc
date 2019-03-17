#include <iostream>
#include <sglib/workspace/WorkSpace.hpp>
int main(int argc, char **argv) {

    WorkSpace ws;
    ws.load_from_disk("nano10x_rg_lmpLR10x_mapped.bsgws");

    ws.long_read_mappers[0].read_filtered_mappings("fm_10K3.bsgfrm");
    ws.long_read_mappers[0].update_indexes();
//    ws.long_read_mappers[0].improve_read_filtered_mappings(100166);

//    for (int i = 35729; i < 38174; i++) {
//        ws.long_read_mappers[0].create_read_path()
//    }

    for (int i = 37174; i < 39174; i++) {
//    ws.long_read_mappers[0].improve_read_filtered_mappings(i, true);
    ws.long_read_mappers[0].create_read_path(i);
    }

    return 0;
}