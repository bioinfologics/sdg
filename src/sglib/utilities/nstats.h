//
// Created by Ben Ward (EI) on 12/04/2018.
//

#ifndef BSG_N_STATS_H
#define BSG_N_STATS_H

#include <numeric>
#include <algorithm>

unsigned long long NX(std::vector<uint32_t>& sizes, const int N) {
    const auto total_size = std::accumulate(sizes.begin(), sizes.end(), 0);
    const unsigned long long limit = ((double)total_size / 100) * N;
    std::sort(sizes.begin(), sizes.end(), std::greater<uint32_t>());
    unsigned long long NX = 0, accum = 0;
    for (const auto& current_size : sizes) {
        NX = current_size;
        accum += current_size;
        if (accum >= limit) break;
    }
    return NX;
}

#endif //BSG_N_STATS_H
