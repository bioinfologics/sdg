//
// Created by Luis Yanes (EI) on 14/08/2018.
//

#ifndef BSG_VERSIONING_HPP
#define BSG_VERSIONING_HPP

#include <cstdint>
typedef uint16_t sdgVersion_t;
typedef uint16_t sdgMagic_t;

static const sdgMagic_t SDG_MAGIC = 0x05D6;
static const sdgVersion_t SDG_VN = 0x0003;

enum SDG_FILETYPE : uint16_t{
    WS_FT,
    KCI_FT,
    PairedDS_FT,
    LinkedDS_FT,
    LongDS_FT,
    PathDS_FT,
    PairedMap_FT,
    LinkedMap_FT,
    LongMap_FT,
    HLAP_FT,
    HLAF_FT,
    NUM_TYPES
};
#endif //BSG_VERSIONING_HPP
