//
// Created by Luis Yanes (EI) on 17/11/2017.
//

#ifndef SG_COMMON_H_H
#define SG_COMMON_H_H
struct ReaderStats {
    uint64_t totalRecords = 0;
    uint64_t filteredRecords = 0;
    uint64_t totalLength = 0;
    uint64_t filteredLength = 0;
};
#endif //SG_COMMON_H_H
