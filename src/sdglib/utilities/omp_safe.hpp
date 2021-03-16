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
#include <algorithm>
namespace sdglib {
    template<typename RItr>
    void inline sort(RItr first, RItr last) {
        std::sort(first, last);
    }
    template<typename RItr, class Comp>
    void inline sort(RItr first, RItr last, Comp comp_fn) {
        std::sort(first, last, comp_fn);
    }
}
//#include <parallel/algorithm>
//namespace sdglib {
//    template<typename RItr>
//    void inline sort(RItr first, RItr last) {
//        __gnu_parallel::sort(first, last);
//    }
//    template<typename RItr, class Comp>
//    void inline sort(RItr first, RItr last, Comp comp_fn) {
//        __gnu_parallel::sort(first, last, comp_fn);
//    }
//
//
//}
#else
constexpr static inline int omp_get_max_threads();
constexpr static inline int omp_get_num_threads();
static inline int omp_get_thread_num();
#endif

#ifndef _OPENMP
constexpr int omp_get_max_threads() {return 1u;}
constexpr int omp_get_num_threads() {return 1u;}
int omp_get_thread_num() {return 0u;}

#include <algorithm>
namespace sdglib {
    template<typename RItr>
    void inline sort(RItr first, RItr last) {
        std::sort(first, last);
    }
    template<typename RItr, class Comp>
    void inline sort(RItr first, RItr last, Comp comp_fn) {
        std::sort(first, last, comp_fn);
    }
}

#endif

#endif //BSG_OMP_SAFE_HPP
