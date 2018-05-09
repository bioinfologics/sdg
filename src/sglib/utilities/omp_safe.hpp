//
// Created by Luis Yanes (EI) on 23/03/2018.
//

#ifndef BSG_OMP_SAFE_HPP
#define BSG_OMP_SAFE_HPP


/**
 * Defines some commonly used OMP header dependent functions in order to have
 * compatibility with compilers that do not support OpenMP
 */
#ifdef _OPENMP
#include <omp.h>
#include <parallel/algorithm>
#else
static inline int omp_get_max_threads();
static inline int omp_get_thread_num();
#endif

#ifndef _OPENMP
int omp_get_max_threads() {return 1u;}
int omp_get_thread_num() {return 0u;}
#endif

#endif //BSG_OMP_SAFE_HPP
