//
// Created by Ben Ward (EI) on 15/03/2018.
//

#ifndef BSG_MATRIX_H
#define BSG_MATRIX_H


#include <vector>
#include <algorithm>
#include <iterator>
#include <utility>

namespace sglib {
    namespace utilities {
        namespace matrices {

            template<class T>
            class matrix {
            public:
                // misc types
                using data_type  = std::vector<T>;
                using value_type = typename std::vector<T>::value_type;
                using size_type  = typename std::vector<T>::size_type;
                // ref
                using reference       = typename std::vector<T>::reference;
                using const_reference = typename std::vector<T>::const_reference;
                // iter
                using iterator       = typename std::vector<T>::iterator;
                using const_iterator = typename std::vector<T>::const_iterator;
                // reverse iter
                using reverse_iterator       = typename std::vector<T>::reverse_iterator;
                using const_reverse_iterator = typename std::vector<T>::const_reverse_iterator;

                // empty construction
                matrix() = default;

                // default-insert rows*cols values
                matrix(size_type rows, size_type cols)
                        : m_rows (rows), m_cols(cols), m_data(rows * cols) {}

                // copy initialized matrix rows*cols
                matrix(size_type rows, size_type cols, const_reference val)
                        : m_rows (rows), m_cols(cols), m_data(rows * cols, val) {}

                matrix(size_type rows, size_type cols, const data_type& data)
                        : m_rows (rows), m_cols(cols), m_data(data) {}

                // 1d-iterators

                iterator begin() { return m_data.begin(); }

                iterator end() { return m_data.end(); }

                const_iterator begin() const { return m_data.begin(); }

                const_iterator end() const { return m_data.end(); }

                const_iterator cbegin() const { return m_data.cbegin(); }

                const_iterator cend() const { return m_data.cend(); }

                reverse_iterator rbegin() { return m_data.rbegin(); }

                reverse_iterator rend() { return m_data.rend(); }

                const_reverse_iterator rbegin() const { return m_data.rbegin(); }

                const_reverse_iterator rend() const { return m_data.rend(); }

                const_reverse_iterator crbegin() const { return m_data.crbegin(); }

                const_reverse_iterator crend() const { return m_data.crend(); }

                // element access (row major indexation)
                reference operator()(size_type const row, size_type const column) {
                    return m_data[m_cols * row + column];
                }

                const_reference operator()(size_type const row, size_type const column) const {
                    return m_data[m_cols * row + column];
                }

                reference at(size_type const row, size_type const column) {
                    return m_data.at(m_cols * row + column);
                }

                const_reference at(size_type const row, size_type const column) const {
                    return m_data.at(m_cols * row + column);
                }

                // resizing
                void resize(size_type new_rows, size_type new_cols) {
                    // new matrix new_rows times new_cols
                    matrix tmp(new_rows, new_cols);
                    // select smaller row and col size
                    auto mc = std::min(m_cols, new_cols);
                    auto mr = std::min(m_rows, new_rows);
                    for (size_type i(0U); i < mr; ++i) {
                        // iterators to begin of rows
                        auto row = begin() + i * m_cols;
                        auto tmp_row = tmp.begin() + i * new_cols;
                        // move mc elements to tmp
                        std::move(row, row + mc, tmp_row);
                    }
                    // move assignment to this
                    *this = std::move(tmp);
                }

                // size and capacity
                size_type size() const { return m_data.size(); }

                size_type max_size() const { return m_data.max_size(); }

                bool empty() const { return m_data.empty(); }

                // dimensionality
                size_type rows() const { return m_rows; }

                size_type cols() const { return m_cols; }

                // data swapping
                void swap(matrix &rhs) {
                    using std::swap;
                    m_data.swap(rhs.m_data);
                    swap(m_rows, rhs.m_rows);
                    swap(m_cols, rhs.m_cols);
                }

                void show_matrix() const {
                    std::cout << m_rows << " x " << m_cols << " matrix:" << std::endl;
                    for(size_type r = 0; r < m_rows; ++r) {
                        for(size_type c = 0; c < m_cols; ++c) {
                            printf("%2d, ", at(r, c));
                        }
                        std::cout << std::endl;
                    }
                }

            private:
                // content
                size_type m_rows{0u};
                size_type m_cols{0u};
                data_type m_data{};
            };

            template<class T>
            void swap(matrix <T> &lhs, matrix <T> &rhs) {
                lhs.swap(rhs);
            }

            template<class T>
            bool operator==(matrix <T> const &a, matrix <T> const &b) {
                if (a.rows() != b.rows() || a.cols() != b.cols()) {
                    return false;
                }
                return std::equal(a.begin(), a.end(), b.begin(), b.end());
            }

            template<class T>
            bool operator!=(matrix <T> const &a, matrix <T> const &b) {
                return !(a == b);
            }
        }
    }
}



#endif //BSG_MATRIX_H
