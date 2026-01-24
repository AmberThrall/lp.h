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

#include <algorithm>
#include <cstdint>
#include <vector>
#include <ostream>
#include <iomanip>
#include <stdexcept>

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

        /// Returns the n x n identity matrix
        static Matrix identity(size_t n) {
            Matrix i(n, n);
            for (size_t k = 0; k < n; ++k) {
                i(k, k) = 1;
            }
            return i;
        }

        /// Returns the augmented matrix [a | b]
        static Matrix augment(const Matrix& a, const Matrix& b) {
            if (a.rows() != b.rows()) {
                throw std::invalid_argument("augment: dimension mismatch");
            }

            Matrix aug(a.rows(), a.cols() + b.cols());
            for (const auto e : a) {
                aug(e.row, e.col) = e.value;
            }

            for (const auto e : b) {
                aug(e.row, e.col + a.cols()) = e.value;
            }

            return aug;
        }

        /// Returns the number of rows
        size_t rows() const { return n_rows; }

        /// Returns the number of columns
        size_t cols() const { return n_cols; }

        /// Returns the number of non-zero entries
        size_t non_zeros() const { return value.size(); }

        bool square() const { return rows() == cols(); }

        /// Returns a submatrix formed by a list of columns
        Matrix submatrix_cols(std::vector<size_t> subcols) const {
            Matrix ret(n_rows, subcols.size());
            for (size_t i = 0; i < value.size(); ++i) {
                for (size_t j = 0; j < subcols.size(); ++j) {
                    if (subcols[j] == col[i]) {
                        ret(row[i], j) = value[i];
                    }
                }
            }
            return ret;
        }

        /// Returns a submatrix formed by a list of rows
        Matrix submatrix_rows(std::vector<size_t> subrows) const {
            Matrix ret(subrows.size(), cols());
            for (size_t i = 0; i < value.size(); ++i) {
                for (size_t j = 0; j < subrows.size(); ++j) {
                    if (subrows[j] == row[i]) {
                        ret(j, col[i]) = value[i];
                    }
                }
            }
            return ret;
        }

        /// Swaps two rows.
        void swap_rows(size_t r1, size_t r2) {
            for (size_t i = 0; i < value.size(); ++i) {
                if (row[i] == r1) { row[i] = r2; }
                else if (row[i] == r2) { row[i] = r1; }
            }
        }

        /// Scales a row by `s`
        void scale_row(size_t r, Number s) {
            for (size_t i = 0; i < value.size(); ++i) {
                if (row[i] == r) { value[i] *= s; } 
            }
        }

        /// Performs basic ERO R1 <- R1 + s*R2
        void add_rows(size_t r1, size_t r2, Number s) {
            for (size_t i = 0; i < value.size(); ++i) {
                if (row[i] == r2) { 
                    (*this)(r1, col[i]) += s * value[i];
                } 
            }
        }

        /// Computes the matrix's row reduced form
        void rref() {
            size_t lead = 0;
            for (size_t r = 0; r < rows(); ++r) {
                if (lead >= cols()) { return; }

                // Find the pivot
                size_t i = r;
                while ((*this)(i, lead) == 0) {
                    i += 1;
                    if (i == rows()) {
                        i = r;
                        lead += 1;
                        if (lead == cols()) { return; }
                    }
                }

                swap_rows(i, r);
                scale_row(r, 1/(*this)(r, lead));
                for (size_t i = 0; i < rows(); ++i) {
                    if (i == r) continue;
                    add_rows(i, r, -(*this)(i, lead));
                }

                 lead += 1;
            }
        }

        /// Computes a matrix's inverse. Warning: may produce incorrect results if the inverse does not exist
        Matrix inverse() {
            if (!square()) {
                throw std::invalid_argument("inverse: matrix is not square");
            }

            std::vector<size_t> c;
            for (size_t i = 0; i < cols(); ++i) { c.push_back(cols() + i); }

            Matrix aug = Matrix::augment(*this, Matrix::identity(rows()));
            aug.rref();
            return aug.submatrix_cols(c);
        }

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
            if (r >= rows() || c >= cols()) {
                throw std::invalid_argument("index out-of-bounds");
            }
            
            for (size_t i = 0; i < value.size(); ++i) {
                if (row[i] == r && col[i] == c) {
                    return value[i];
                }
            }
            return 0;
        }

        /// Access a reference to entry at (r,c), inserts a 0 if entry does not exist
        Number& operator()(std::size_t r, std::size_t c) {
            if (r >= rows() || c >= cols()) {
                throw std::invalid_argument("index out-of-bounds");
            }

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

        /// Matrix multiplication
        Matrix operator*(const Matrix& rhs) const {
            if (cols() != rhs.rows()) {
                throw std::invalid_argument("multiply: dimension mismatch");
            }

            Matrix C(rows(), rhs.cols());            

            // Group rhs by rows
            std::vector<std::vector<size_t>> brow(rhs.rows());
            for (size_t k = 0; k < rhs.value.size(); ++k) {
                brow[rhs.row[k]].push_back(k);
            }

            for (size_t a = 0; a < value.size(); ++a) {
                size_t r = row[a];
                size_t k = col[a];

                for (size_t b : brow[k]) {
                    size_t j = rhs.col[b];
                    C(r, j) += value[a] * rhs.value[b];
                }
            }

            return C;
        }

        void operator+=(const Matrix& rhs) {
            if (cols() != rhs.cols() || rows() != rhs.rows()) {
                throw std::invalid_argument("add: dimension mismatch");
            }

            // Group rhs by rows
            for (const auto& e : rhs) {
                (*this)(e.row, e.col) += e.value; 
            }
        }

        void operator-=(const Matrix& rhs) {
            if (cols() != rhs.cols() || rows() != rhs.rows()) {
                throw std::invalid_argument("subtract: dimension mismatch");
            }

            // Group rhs by rows
            for (const auto& e : rhs) {
                (*this)(e.row, e.col) -= e.value; 
            }
        }

        Matrix operator*(const Number& rhs) const {
            Matrix C(*this);
            for (size_t i = 0; i < C.value.size(); ++i) {
                C.value[i] *= rhs;
            }
            return C;
        }
    private:
        std::size_t n_rows, n_cols;
        std::vector<Number> value;
        std::vector<size_t> col;
        std::vector<size_t> row;
    };

    inline std::ostream& operator<<(std::ostream& os, const Matrix& m) {
        size_t min_w = 8;

        for (size_t r = 0; r < m.rows(); ++r) {
            for (size_t c = 0; c < m.cols(); ++c) {
                os << std::setw(min_w) << std::fixed << std::setprecision(2) << m(r, c);
            }
            os << std::endl;
        }

        return os;
    }

    enum class SolutionStatus { kOptimal, kUnbounded, kInfeasible, kFeasible, kAborted };


    /// Performs revised simplex method to solve the LP:
    /// min c^Tx
    /// s.t. Ax <= b
    ///       x >= 0
    class RevisedSimplex {
    public:
        struct State {
            Matrix& A;
            Matrix& b;
            Matrix& c;
            Matrix& x;
            std::vector<size_t> bv;
            std::vector<size_t> nbv;
            Matrix Binv;
            SolutionStatus status;
        };
        
        struct Solution {
            std::vector<Number> x;
            Number z;
            SolutionStatus status;
        };

        typedef void (*HookFn)(State& state);

        RevisedSimplex(Matrix& a, Matrix& b, Matrix& c) : A(a), b(b), c(c), hook(nullptr) {
            A = Matrix::augment(a, Matrix::identity(a.rows()));
            c = Matrix::augment(c, Matrix(1,a.rows()));
        }

        /// Sets a hook function that is called each iteration
        void set_hook(HookFn fn) { hook = fn; }

        Solution solve() {
            Matrix c_orig(c);

            // Always use the slack variables as starting basis variables
            std::vector<size_t> bv;
            std::vector<size_t> nbv;
            for (size_t i = 0; i < A.cols(); ++i) {
                if (i >= A.cols() - A.rows()) { bv.push_back(i); } 
                else { nbv.push_back(i); }
            }
            Matrix Binv = Matrix::identity(bv.size());
            Matrix x(A.cols(), 1);
            for (const auto& e: b) { x(e.row + nbv.size(), e.col) = e.value; }

            State state = { A, b, c, x, bv, nbv, Binv, SolutionStatus::kFeasible };

            while (state.status == SolutionStatus::kFeasible) {
                step(state);
            }

            // Create the solution
            Number z = (c_orig * state.x)(0,0);
            std::vector<Number> x_soln(A.cols());
            for (const auto& e: state.x) { x_soln[e.row] = e.value; }

            return Solution { x_soln, z, state.status };
        }
    private:
        void step(State& state) {
            // Check if optimal
            Matrix cN = state.c.submatrix_cols(state.nbv);
            Matrix cB = state.c.submatrix_cols(state.bv);
            cN -= cB * state.Binv * state.A.submatrix_cols(state.nbv);

            bool optimal = true;
            for (const auto& e : cN) { 
                if (e.value < 0) { optimal = false; break; }
            }
            if (optimal) { state.status = SolutionStatus::kOptimal; }

            // Call the hook
            if (hook != nullptr) {
                hook(state);
            }
            if (state.status != SolutionStatus::kFeasible) { return; }
                        
            // Determine the entering variable
            size_t enter_var = 1;
            Number best = 0;
            for (const auto& e : cN) {
                if (e.value < best) {
                    enter_var = state.nbv[e.col];
                    best = e.value;
                }
            }
 
            // Check if LP is unbounded
            Matrix dB = state.Binv * state.A.submatrix_cols({enter_var}) * -1;
            bool unbounded = true;
            for (const auto& e : dB) {
                if (e.value < 0) { unbounded = false; break; }
            }

            if (unbounded) {
                state.status = SolutionStatus::kUnbounded;
                return;
            }

            // Perform minimum ratio test for leaving variable
            Matrix xB = state.x.submatrix_rows(state.bv); 
            
            Number min_ratio = 0;
            size_t min_ratio_winner = 0;
            for (size_t i = 0; i < state.bv.size(); ++i) {
                Number x = xB(i, 0);
                Number d = dB(i, 0);
                if (x > 0) {
                    Number ratio = -x / d;
                    if (ratio < min_ratio || min_ratio == 0) {
                        min_ratio = ratio;
                        min_ratio_winner = i;
                    }
                }
            }

            if (min_ratio == 0) {
                state.status = SolutionStatus::kUnbounded;
                return;
            }       

            Matrix d(state.x.rows(), 1);
            for (const auto& e: dB) { d(state.bv[e.row], 0) = e.value; }
            d(enter_var, 0) = 1;
            state.x += d * min_ratio;
            
            size_t leaving_var = state.bv[min_ratio_winner];

            state.bv[min_ratio_winner] = enter_var;
            std::sort(state.bv.begin(), state.bv.end());

            for (size_t i = 0; i < state.nbv.size(); ++i) {
                if (state.nbv[i] == enter_var) { state.nbv[i] = leaving_var; break; }
            }
            std::sort(state.nbv.begin(), state.nbv.end());

            // Compute the new Binv
            Matrix aug = Matrix::augment(state.Binv, dB * -1);

            // We want last column of aug to be e_j, where j = enter_var
            aug.scale_row(enter_var, -1/dB(enter_var, 0));
            for (const auto& e : dB) {
                if (e.row == enter_var) { continue; }
                aug.add_rows(e.row, enter_var, e.value);
            }

            std::vector<size_t> subcols;
            for (size_t i = 0; i < state.Binv.cols(); ++i) { subcols.push_back(i); }
            state.Binv = aug.submatrix_cols(subcols);
        }

        Matrix& A;
        Matrix& b;
        Matrix& c;
        HookFn hook;
    };
}
