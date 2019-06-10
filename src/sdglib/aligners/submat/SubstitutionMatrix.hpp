//
// Created by Ben Ward (EI) on 15/03/2018.
//

#ifndef BSG_SUBSTITUTIONMATRIX_H
#define BSG_SUBSTITUTIONMATRIX_H

#include <iostream>
#include <sdglib/utilities/matrix.hpp>

namespace sdglib {
    namespace alignment {
        namespace scoring {

            // TODO: Make substitution matrix work for AA sequences in future too, not just DNA.

            template<class T>
            class SubstitutionMatrix : public utilities::matrices::matrix<T> {
            public:

                using size_type  = typename utilities::matrices::matrix<T>::size_type;
                using data_type  = typename utilities::matrices::matrix<T>::data_type;
                using value_type = typename utilities::matrices::matrix<T>::value_type;

                SubstitutionMatrix(size_type rows, size_type cols, const data_type& data)
                        : utilities::matrices::matrix<T>(rows, cols, data) {}

                value_type at(char row, char column) const {
                    auto X = (size_t) std::toupper(row) - 65;
                    auto Y = (size_t) std::toupper(column) - 65;
                    size_type r = DNA_CHAR_TO_SIZE_T[X];
                    size_type c = DNA_CHAR_TO_SIZE_T[Y];
                    return utilities::matrices::matrix<T>::at(r, c);
                }

            private:
                static constexpr size_type DNA_CHAR_TO_SIZE_T[] = {  0, 13,  1, 12, 15,
                                                                    15,  3, 10, 15, 15,
                                                                    11, 15,  2, 14, 15,
                                                                    15, 15,  4,  5,  7,
                                                                    15,  6,  8, 15,  9 };
            };

            template<class T> constexpr typename SubstitutionMatrix<T>::size_type SubstitutionMatrix<T>::DNA_CHAR_TO_SIZE_T[];

            const SubstitutionMatrix<int> EDNAFULL {15, 15, { 5, -4,  1, -4,  1, -4, -1, -4,  1, -4, -1, -4, -1, -4, -2,
                                                             -4,  5,  1, -4, -4,  1, -1, -4, -4,  1, -1, -4, -4, -1, -2,
                                                              1,  1, -1, -4, -2, -2, -1, -4, -2, -2, -1,  4, -3, -3, -1,
                                                             -4, -4, -4,  5,  1,  1, -1, -4, -4, -4, -4,  1, -1, -1, -2,
                                                              1, -4, -2,  1, -1, -2, -1, -4, -2, -4, -3, -2, -1, -3, -1,
                                                             -4,  1, -2,  1, -2, -1, -1, -4, -4, -2, -3, -2, -3, -1, -1,
                                                             -1, -1, -1, -1, -1, -1, -1, -4, -3, -3, -2, -3, -2, -2, -1,
                                                             -4, -4, -4, -4, -4, -4, -4,  5,  1,  1, -1,  1, -1, -1, -2,
                                                              1, -4, -2, -4, -2, -4, -3,  1, -1, -2, -1, -2, -1, -3, -1,
                                                             -4,  1, -2, -4, -4, -2, -3,  1, -2, -1, -1, -2, -3, -1, -1,
                                                             -1, -1, -1, -4, -3, -3, -2, -1, -1, -1, -1, -3, -2, -2, -1,
                                                             -4, -4, -4,  1, -2, -2, -3,  1, -2, -2, -3, -1, -1, -1, -1,
                                                             -1, -4, -3, -1, -1, -3, -2, -1, -1, -3, -2, -1, -1, -2, -1,
                                                             -4, -1, -3, -1, -3, -1, -2, -1, -3, -1, -2, -1, -2, -1, -1,
                                                             -2, -2, -1, -2, -1, -1, -1, -2, -1, -1, -1, -1, -1, -1, -1 }};


        }
    }
}


#endif //BSG_SUBSTITUTIONMATRIX_H
