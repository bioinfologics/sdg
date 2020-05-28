//
// Created by Ben J. Ward (EI) on 27/05/2020.
//

#ifndef SDG_KMERITERATOR_H
#define SDG_KMERITERATOR_H

#pragma once

// QUESTION: Should SDG get it's own sequence classes that aren't just strings?
// Check BioSequences.jl for example of what I can do.

/*
class KmerIterator : public std::iterator<std::input_iterator_tag, uint64_t> {
    // C++11 (explicit aliases)
    using iterator_category = std::input_iterator_tag;
    using value_type = uint64_t;
    using reference = value_type const&;
    using pointer = value_type const*;
    using difference_type = ptrdiff_t;

    // C++11
    explicit operator bool() const { return !done; }
    reference operator*() const { return fwmer; }
    pointer operator->() const { return &fwmer; }

    KmerIterator& operator++() {
        assert(!done);
        char nt = seq[pos];
        value_type fbits = kmerbits[nt]; // Get bits for the nucleotide.
        value_type rbits = ~fbits & (value_type) 3;
        fwmer = (fwkmer << 2 | fbits);
        rvmer = (rvkmer >> 2) | (rbits << ((k - 1) * 2));
        pos == stop && done = true;
        pos++;
        return *this;
    }

    KmerIterator(const std::string& seq, const uint8_t k) : seq(seq), k(k), done(false), pos(0), fwmer(0), rvmer(0) {
        size_t filled = 0;
        size_t i = 0;
        while(i < seq.length()){
            nt = reinterpret(UInt8, inbounds_getindex(it.seq, i))
            @inbounds fbits = UT(kmerbits[nt + 1])
            rbits = ~fbits & typeof(fbits)(0x03)
            fwkmer = (fwkmer << 0x02 | fbits)
            rvkmer = (rvkmer >> 0x02) | (UT(rbits) << unsigned(offset(T, 1)))
            filled += 1
            if filled == ksize(T)
            return MerIterResult(1, T(fwkmer), T(rvkmer)), (i, fwkmer, rvkmer)
            end
            i += 1
        }
    }
private:
    const std::string& seq;
    const uint8_t k;
    bool done;
    size_t pos;
    value_type fwmer;
    value_type rvmer;
};
*/

// I want this function to be run at compile time to create the constant static
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
template<int K> // I'm not sure if K should be a type parameter - it is for my Kmer types in BioJulia.
// It means some stuff and constants get worked out at compile time e.g. see value_type below,
// and it stops people mixing kmers of different K and shit like that.
// Not sure if as valid for C++ as C++ is not JIT and much less interactive.
class Kmers {
    // C++11 (explicit aliases for input iterator).
    using iterator_category = std::input_iterator_tag;
    using value_type = typename std::conditional<(K > 32), uint64_t, __uint128_t>::type;
    using reference = value_type const&;
    using pointer = value_type const*;
    using difference_type = ptrdiff_t;

    // Static array for the conversion from letter to character.
    static const std::array<value_type, 256> nuc_to_bin = make_nuc2bin_tbl<value_type>();
public:
    Kmers(const std::string& seq) : sequence(seq) {
        // TODO: Add some validation of the string?
    }

    class iterator: public std::iterator<iterator_category, value_type, difference_type, pointer, reference> {
    public:
        explicit iterator(const Kmers& parent, const size_t kmerend) : parent(parent), pos(0) {
            stop = parent.sequence.length();
            for(int i = 0; i < stop; ++i) {
                this++;
            }
        }
        reference operator*() const { return std::min(fwmer, rvmer); }
        pointer operator->() const { return std::min(fwmer, rvmer); }
        iterator& operator++() {
            // Check we are not done.
            // Then add the base pointed to by nextbase to rvmer and fwmer.
            // Then increment nextbase
            assert(!done);
            char nucleotide = parent.sequence[nextbase];
            // Get bits for the nucleotide.
            // The array should be a static and constant array.
            value_type fbits = parent.nuc_to_bin[nucleotide];
            value_type rbits = ~fbits & (value_type) 3;
            fwmer = (fwmer << 2 | fbits);
            rvmer = (rvmer >> 2) | (rbits << ((K - 1) * 2));
            nextbase == lastbase && done = true;
            nextbase++;
            return *this;
        }
    private:
        const Kmers& parent;
        size_t nextbase;
        size_t lastbase;
        value_type fwmer;
        value_type rvmer;
        bool done;
    };
    iterator begin() { return iterator(this, K); }
    iterator end() {}
private:
    std::string& sequence;
};

#endif //SDG_KMERITERATOR_H
