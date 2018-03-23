//
// Created by Luis Yanes (EI) on 23/03/2018.
//

#ifndef BSG_OMP_SAFE_HPP
#define BSG_OMP_SAFE_HPP

#ifndef __APPLE__
#include <omp.h>
#include <parallel/algorithm>
#else
static inline int omp_get_max_threads();
static inline int omp_get_thread_num();
#endif

#ifdef __APPLE__
int omp_get_max_threads() {return 1u;}
int omp_get_thread_num() {return 0u;}
#endif

#endif //BSG_OMP_SAFE_HPP
