//
// Created by Luis Yanes (EI) on 14/08/2018.
//

#ifndef BSG_VERSIONING_HPP
#define BSG_VERSIONING_HPP

typedef uint16_t bsgVersion_t;
typedef uint16_t bsgMagic_t;

static const bsgMagic_t BSG_MAGIC = 0x0B56;
static const bsgVersion_t BSG_VN = 0x0001;

enum BSG_FILETYPE : uint16_t{
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
