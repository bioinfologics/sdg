//
// Created by Ben Ward (EI) on 15/03/2018.
//

#ifndef BSG_SMITHWATERMAN_H
#define BSG_SMITHWATERMAN_H

#include <numeric>
#include <string>
#include <vector>
#include <sdglib/utilities/matrix.hpp>
#include <sdglib/aligners/submat/SubstitutionMatrix.hpp>

namespace sdglib {
    namespace alignment {
        namespace algorithms {

            typedef uint8_t trace_t;

            const trace_t TRACE_NONE = 0b00000;
            const trace_t TRACE_MATCH = 0b00001;
            const trace_t TRACE_DELETE = 0b00010;
            const trace_t TRACE_INSERT = 0b00100;
            const trace_t TRACE_EXTDEL = 0b01000;
            const trace_t TRACE_EXTINS = 0b10000;

            using namespace sdglib::utilities::matrices;
            using namespace sdglib::alignment::scoring;


            template<class T>
            class SmithWaterman {
            public:
                using size_type = matrix<trace_t>::size_type;

                SmithWaterman(size_type m, size_type n) : traces({m + 1, n + 1, 0xff}), H(m + 1), E(m) { }

                void show_trace() const { traces.show_matrix(); }
                void show_H() const { for (const auto& i : H) std::cout << i << ", "; std::cout << std::endl; }
                void show_E() const { for (const auto& i : E) std::cout << i << ", "; std::cout << std::endl; }

                std::tuple<T, std::tuple<int, int>> run(const std::string &seqA, const std::string &seqB,
                                                        T gap_open, T gap_extend,
                                                        const SubstitutionMatrix<T> &submat = EDNAFULL) {
                    return run(seqA, seqB, gap_open, gap_extend, gap_open, gap_extend, submat);
                };

                std::tuple<T, std::tuple<int, int>> run(const std::string &seqA, const std::string &seqB,
                                                        T gap_open_a, T gap_extend_a,
                                                        T gap_open_b, T gap_extend_b,
                                                        const SubstitutionMatrix<T> &submat = EDNAFULL) {

                    std::string::size_type lengthSeqA = seqA.length();
                    std::string::size_type lengthSeqB = seqB.length();

                    ensure_room(lengthSeqA, lengthSeqB);

                    auto gap_init_a = gap_open_a + gap_extend_a;
                    auto gap_init_b = gap_open_b + gap_extend_b;

                    H[0] = (T) 0;
                    traces(0, 0) = TRACE_NONE;
                    for (size_t i = 0; i < lengthSeqA; ++i) {
                        H[i + 1] = 0;
                        E[i] = H[i + 1] + gap_init_a;
                        traces(i + 1, 0) = TRACE_NONE;
                        if (lengthSeqB >= 1) {
                            traces(i + 1, 1) = TRACE_NONE;
                        }
                    }

                    T best_score = H[0];
                    std::tuple<int, int> best_endpos{0, 0};

                    for (size_t j = 0; j < lengthSeqB; ++j) {
                        char b_j = seqB[j];
                        T h_diag = H[0];
                        H[0] = 0;
                        T f = H[0] + gap_init_b;
                        trace_t ft = TRACE_NONE;
                        traces(0, j + 1) = TRACE_NONE;

                        for (size_t i = 0; i < lengthSeqA; ++i) {
                            T e = E[i];
                            T g = h_diag + submat.at(b_j, seqA[i]);
                            T h = std::max(std::initializer_list<T>({(T) 0, e, f, g}));
                            h_diag = H[i + 1];
                            H[i + 1] = h;
                            trace_t t = traces(i + 1, j + 1) | ft;
                            if (e == h) t |= TRACE_DELETE;
                            if (f == h) t |= TRACE_INSERT;
                            if (g == h) t |= TRACE_MATCH;
                            traces(i + 1, j + 1) = t;
                            if (h >= best_score) {
                                best_score = h;
                                best_endpos = {i, j};
                            }
                            // next E
                            if (j != (lengthSeqB - 1)) {
                                T e_prime = e + gap_extend_a;
                                e = std::max(e_prime, h + gap_init_a);
                                E[i] = e;
                                trace_t et = TRACE_NONE;
                                if (e == e_prime) et |= TRACE_EXTDEL;
                                traces(i + 1, j + 2) = et;
                            }
                            // next F
                            T f_prime = f + gap_extend_b;
                            f = std::max(f_prime, h + gap_init_b);
                            ft = TRACE_NONE;
                            if (f == f_prime) ft |= TRACE_EXTINS;
                        }
                    }

                    return std::tuple<T, std::tuple<int, int>> {best_score, best_endpos};
                }

            private:
                matrix<trace_t> traces;
                std::vector<T> H, E;

                void ensure_room(size_type m, size_type n) {
                    if (traces.size() != ((m + 1) * (n + 1))) {
                        traces.resize(m + 1, n + 1);
                        H.resize(m + 1);
                        E.resize(m);
                    }
                }
            };
        }
    }
}
#endif //BSG_SMITHWATERMAN_H
