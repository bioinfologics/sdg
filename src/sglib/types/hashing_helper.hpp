//
// Created by Luis Yanes (EI) on 26/03/2018.
//

#ifndef BSG_HASHING_HELPER_HPP
#define BSG_HASHING_HELPER_HPP
#include <tuple>

/**
 * Takes each of the elements in a tuple and combines them in a single hash value of size_t length
 */
namespace sglib {
    namespace
    {
        template <typename TT>
        struct hash
        {
            size_t
            operator()(TT const& tt) const
            {
                return std::hash<TT>()(tt);
            }
        };


        template <class T>
        inline void hash_combine(std::size_t& seed, T const& v)
        {
            seed ^= sglib::hash<T>()(v) + 0x9e3779b9 + (seed<<6) + (seed>>2);
        }

        template <class TP,
                size_t index = std::tuple_size<TP>::value - 1>
        struct hashvalueimpl
        {
            static void apply(size_t& seed, TP const& tuple)
            {
                hashvalueimpl<TP, index-1>::apply(seed, tuple);
                hash_combine(seed, std::get<index>(tuple));
            }
        };

        template <class TP>
        struct hashvalueimpl<TP,0>
        {
            static void apply(size_t& seed, TP const& tuple)
            {
                hash_combine(seed, std::get<0>(tuple));
            }
        };
    }

    template <class ... TT>
    struct hash<std::tuple<TT...>>
    {
        size_t
        operator()(std::tuple<TT...> const& tuple) const
        {
            size_t seed = 0;
            hashvalueimpl<std::tuple<TT...> >::apply(seed, tuple);
            return seed;
        }

    };
}
#endif //BSG_HASHING_HELPER_HPP
