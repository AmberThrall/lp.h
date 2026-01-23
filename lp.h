// MIT License
// 
// Copyright (c) 2026 Amber R. Thrall
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in all
// copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.


/**
 *
 */

#pragma once

#define LP_H_VERSION_MAJOR "0"
#define LP_H_VERSION_MINOR "1"
#define LP_H_VERSION_PATCH "0"
#define LP_H_VERSION "v" LP_H_VERSION_MAJOR "." LP_H_VERSION_MINOR "." LP_H_VERSION_PATCH

#include <cstdint>
#include <vector>
#include <ostream>
#include <iomanip>

namespace lp {
#ifdef LP_H_USE_FLOAT
    using Number = float;
#else
    using Number = double;
#endif

    class Matrix {
    public:
        /// Creates a m x n all-zeros matrix 
        Matrix(size_t m, size_t n) : n_rows(m), n_cols(n) {}

        /// Returns the number of rows
        size_t rows() const { return n_rows; }

        /// Returns the number of columns
        size_t cols() const { return n_cols; }

        /// Returns the number of non-zero entries
        size_t non_zeros() const { return value.size(); }

        struct Entry {
            size_t row, col;
            Number& value;
        };
        struct ConstEntry {
            size_t row, col;
            const Number& value;
        };

        class iterator {
        public:
            iterator(Matrix* m, size_t idx) : m(m), idx(idx) {}
            Entry operator*() const { return { m->row[idx], m->col[idx], m->value[idx] }; }
            iterator& operator++() { ++idx; return *this; }
            bool operator!=(const iterator& other) const { return idx != other.idx; }
        private:
            Matrix* m;
            size_t idx;
        };

        class const_iterator {
        public:
            const_iterator(const Matrix* m, size_t idx) : m(m), idx(idx) {}
            ConstEntry operator*() const { return { m->row[idx], m->col[idx], m->value[idx] }; }
            const_iterator& operator++() { ++idx; return *this; }
            bool operator!=(const const_iterator& other) const { return idx != other.idx; }
        private:
            const Matrix* m;
            size_t idx;
        };

        iterator begin() { return iterator(this, 0); }
        const_iterator begin() const { return const_iterator(this, 0); }
        iterator end() { return iterator(this, value.size()); }
        const_iterator end() const { return const_iterator(this, value.size()); }

        /// Access the entry at (r,c)
        Number operator()(std::size_t r, std::size_t c) const {
            for (size_t i = 0; i < value.size(); ++i) {
                if (row[i] == r && col[i] == c) {
                    return value[i];
                }
            }
            return 0;
        }

        /// Access a reference to entry at (r,c), inserts a 0 if entry does not exist
        Number& operator()(std::size_t r, std::size_t c) {
            for (size_t i = 0; i < value.size(); ++i) {
                if (row[i] == r && col[i] == c) {
                    return value[i];
                }
            }

            // Insert new entry
            row.push_back(r);
            col.push_back(c);
            value.push_back(0);
            return value.back();
        }
    private:
        std::size_t n_rows, n_cols;
        std::vector<Number> value;
        std::vector<size_t> col;
        std::vector<size_t> row;
    };

    inline std::ostream& operator<<(std::ostream& os, const Matrix& m) {
        size_t min_w = 6;

        for (size_t r = 0; r < m.rows(); ++r) {
            for (size_t c = 0; c < m.cols(); ++c) {
                os << std::setw(min_w) << m(r, c);
            }
            os << std::endl;
        }

        return os;
    }
}
