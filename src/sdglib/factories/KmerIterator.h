//
// Created by Ben J. Ward (EI) on 27/05/2020.
//

#ifndef SDG_KMERITERATOR_H
#define SDG_KMERITERATOR_H

#pragma once

class KmerIterator : public std::iterator<std::input_iterator_tag, uint64_t> {
    // C++11 (explicit aliases)
    using iterator_category = std::input_iterator_tag;
    using value_type = uint64_t;
    using reference = value_type const&;
    using pointer = value_type const*;
    using difference_type = ptrdiff_t;

    // C++11
    explicit operator bool() const { return !done; }
    reference operator*() const { return ij; }
    pointer operator->() const { return &ij; }

private:
    std::string& seq;
    bool done;
    uint64_t fwmer;
};


#endif //SDG_KMERITERATOR_H
