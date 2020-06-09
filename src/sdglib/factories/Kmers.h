//
// Created by Ben J. Ward (EI) on 27/05/2020.
//

#pragma once

#include <array>
#include <limits>
#include <assert.h>
#include <iostream>
#include <cstring>

// I want this function to be able to run at compile time to create the constant static
// arrays that the Kmers class uses to convert chars to binary.
template<typename U>
constexpr std::array<U, 256> make_nuc2bin_tbl(){
    std::array<U, 256> tbl = {4};
    tbl['a'] = tbl['A'] = 0;
    tbl['c'] = tbl['C'] = 1;
    tbl['g'] = tbl['G'] = 2;
    tbl['t'] = tbl['T'] = 3;
    return tbl;
}

// A type representing a "lazy" - if you want to call it that, container of all the k-mers in a sequence.
// Has an input iterator in its scope for iterating from the first k-mer in a sequence to the last.
template <typename U> class Kmers {
protected:
    // C++11 (explicit aliases for input iterator).
    using iterator_category = std::input_iterator_tag;
    using value_type = U;
    using reference = value_type const&;
    using pointer = value_type const*;
    using difference_type = ptrdiff_t;

public:
    // Unwrap std::string and pass to the cstring constructor.
    Kmers(const std::string& seq, const int K) : Kmers(seq.data(), K) {} // Requires C++11.

    Kmers(const char * seq, const int K) : sequence(seq), K(K) {
        // Here I check someone hasn't requested silly K size, before
        // making the local translation table for this container.
        // static asserts below?
        assert(K > 0);
        assert(!std::numeric_limits<value_type>::is_signed);
        // Forgive heavy bracket use - I need to get used to if C++ has the operator precedence I'm used to or not.
        assert(K <= ((std::numeric_limits<value_type>::digits / 2))); // -1 here to keep space for the sign bits.
        translation_table = make_nuc2bin_tbl<value_type>();
    }

    class iterator: public std::iterator<iterator_category, value_type, difference_type, pointer, reference> {
    public:
        explicit iterator(Kmers& p, const size_t offset = 0) : parent(p) {
            position = p.sequence;
            position += offset;
            valid_nucleotides = 0;
            fwmer = 0;
            // TODO: Perhaps move into parent so it is not computed for every iterator or everytime end() is called?
            mask = (value_type(1) << (2 * parent.K)) - 1;
            incorporate_char(position);
        }

        bool has_reached_end() { return *position == '\0'; }

        reference operator*() const {
            return fwmer;
        }

        // Advance the iterator to the next valid kmer in the sequence, or to the end of the string
        // (\0), whichever comes first.
        iterator& operator++() {
            scan_right();
            return *this;
        }
        bool operator<(Kmers::iterator other) {
            // true if both iterators refer to same sequence, and this iterator is at a lower position.
            return parent.sequence == other.parent.sequence && position < other.position;
        }
        bool operator<=(Kmers::iterator other) {
            // true if both iterators refer to same sequence, and this iterator is at a lower position.
            return parent.sequence == other.parent.sequence && position <= other.position;
        }
        bool operator!=(Kmers::iterator other) {
            // true if both iterators refer to the same sequence, and this iterator is not at the same positon.
            return parent.sequence == other.parent.sequence && position != other.position;
        }
        bool operator==(Kmers::iterator other) {
            // true if both iterators refer to the same sequence, and this iterator is not at the same positon.
            return parent.sequence == other.parent.sequence && position == other.position;
        }
    protected:
        const Kmers& parent;
        const char* position;
        value_type fwmer;
        value_type mask;
        size_t valid_nucleotides;

        void scan_right() {
            while(*position != '\0') {
                ++position;
                const char nucleotide = *position;
                // Get bits for the nucleotide.
                // The array should be a static and constant array.
                incorporate_char(nucleotide);
                if(valid_nucleotides >= parent.K) break;
            }
        }

        void incorporate_char(const char c) {
            value_type fbits = parent.translation_table[c];
            value_type rbits = ~fbits & (value_type) 3;
            valid_nucleotides = fbits == 4 ? 0 : valid_nucleotides + 1;
            fwmer = (fwmer << 2 | fbits) & mask;
            //rvmer = ((rvmer >> 2) | (rbits << ((parent.K - 1) * 2))) & mask;
        }
    };
    // An iterator at the first valid Kmer in the sequence.
    iterator begin() {
        iterator it(*this);
        ++it;
        return it;
    }
    // An iterator past all valid Kmers in the sequence (has hit \0 in input string).
    iterator end() {
        iterator it(*this, strlen(sequence));
        return it;
    }
private:
    // std::string& sequence;
    // Doesn't this make it less safe than a std::string reference though?
    // Someone could just free the memory out from underneath us.
    const char* sequence;
    const int K;
    std::array<U, 256> translation_table;
};




// LETS NOT WORRY ABOUT THIS FOR NOW.
template<typename U>
class CanonicalKmers : public Kmers<U> {
    using value_type = typename Kmers<U>::value_type;
    using reference = typename Kmers<U>::reference;
    using pointer = typename Kmers<U>::pointer;
    using Kmers<U>::Kmers; // Inherit constructor.

public:
    class iterator : public Kmers<U>::iterator {
    public:
        explicit iterator(CanonicalKmers& p) : Kmers<U>::iterator(p) {}
        using Kmers<U>::iterator::operator++;
        reference operator*() const {
            // I want to specialize this on template parameter F but am unsure on how to do it properly.
            return std::min(Kmers<U>::iterator::fwmer & Kmers<U>::iterator::mask, Kmers<U>::iterator::rvmer & Kmers<U>::iterator::mask);
        }
        pointer operator->() const {
            // I want to specialize this on template parameter F but am unsure on how to do it properly.
            return std::min(Kmers<U>::iterator::fwmer & Kmers<U>::iterator::mask, Kmers<U>::iterator::rvmer & Kmers<U>::iterator::mask);
        }
    };

    iterator begin() {
        iterator it(*this, 0);
        ++it;
        return it;
    }
    // An iterator past all valid Kmers in the sequence (has hit \0 in input string).
    iterator end() {
        iterator it(*this, 0);
        it.jump_to_end();
        return it;
    }
};
