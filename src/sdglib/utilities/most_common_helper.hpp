//
// Created by Luis Yanes (EI) on 2019-03-18.
//

#ifndef BSG_MOST_COMMON_HELPER_HPP
#define BSG_MOST_COMMON_HELPER_HPP

template<typename A, typename B>
std::pair<B,A> flip_pair(const std::pair<A,B> &p)
{
    return std::pair<B,A>(p.second, p.first);
}

template<typename A, typename B>
std::multimap<B,A> flip_map(const std::map<A,B> &src)
{
    std::multimap<B,A> dst;
    std::transform(src.begin(), src.end(), std::inserter(dst, dst.begin()),
                   flip_pair<A,B>);
    return dst;
}

#endif //BSG_MOST_COMMON_HELPER_HPP
