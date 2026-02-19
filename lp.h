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

#include <cstddef>
#define LP_H_VERSION_MAJOR "0"
#define LP_H_VERSION_MINOR "1"
#define LP_H_VERSION_PATCH "0"
#define LP_H_VERSION "v" LP_H_VERSION_MAJOR "." LP_H_VERSION_MINOR "." LP_H_VERSION_PATCH

#ifdef LP_H_DEBUG
#include <iostream>
#endif
#include <algorithm>
#include <cstdint>
#include <vector>
#include <ostream>
#include <iomanip>
#include <stdexcept>
#include <limits>

namespace lp {
#ifdef LP_H_USE_FLOAT
    using Number = float;
#else
    using Number = double;
#endif
    constexpr Number Infinity = std::numeric_limits<Number>::infinity();
    constexpr Number Eps = 1e-9;

    /// Simple dense vector class
    class Vector {
    public:
        /// Creates a 0-vector
        Vector() : Vector(0) {}

        /// Creates a n-vector filled with v
        Vector(size_t n, Number v=0) {
            value = std::vector<Number>(n, v);
        }

        /// Returns the size of the vector
        size_t size() const { return value.size(); }

        /// Computes the dot product of two vectors
        Number dot(const Vector& other) const {
            if (size() != other.size()) { throw std::invalid_argument("dot product: size mismatch."); }

            Number v = 0; 
            for (size_t i = 0; i < size(); ++i) {
                v += value[i] * other[i];
            }
            return v;
        }

        /// Resizes the vector
        void resize(size_t new_size) {
            value.resize(new_size);
        }

        /// Returns a subvector
        Vector subvector(std::vector<size_t> indices) const {
            Vector ret(indices.size());
            for (size_t i = 0; i < indices.size(); ++i) {
                size_t idx = indices[i];
                if (idx >= size()) { throw std::out_of_range("Vector: out-of-bounds."); }
                ret[i] = value[idx];
            }
            return ret;
        }

        Number& operator[](size_t idx) {
            if (idx >= size()) { throw std::out_of_range("Vector: out-of-bounds."); }
            return value[idx];
        }

        const Number& operator[](size_t idx) const {
            if (idx >= size()) { throw std::out_of_range("Vector: out-of-bounds."); }
            return value[idx];
        }

        void operator+=(const Vector& rhs) {
            if (size() != rhs.size()) { throw std::invalid_argument("vector add: size mismatch"); }

            for (size_t i = 0; i < size(); ++i) {
                value[i] += rhs[i];
            }
        }

        friend Vector operator+(Vector lhs, const Vector& rhs) {
            lhs += rhs;
            return lhs;
        }

        void operator-=(const Vector& rhs) {
            if (size() != rhs.size()) { throw std::invalid_argument("vector add: size mismatch"); }

            for (size_t i = 0; i < size(); ++i) {
                value[i] -= rhs[i];
            }
        }

        friend Vector operator-(Vector lhs, const Vector& rhs) {
            lhs -= rhs;
            return lhs;
        }

        void operator*=(Number rhs) {
            for (size_t i = 0; i < size(); ++i) {
                value[i] *= rhs;
            }
        }

        friend Vector operator*(Vector lhs, Number rhs) {
            lhs *= rhs;
            return lhs;
        }
        friend Vector operator*(Number lhs, Vector rhs) { return rhs * lhs; }
    private:
        std::vector<Number> value;
    };

    inline std::ostream& operator<<(std::ostream& os, const Vector& v) {
        os << "[";
        for (size_t i = 0; i < v.size(); ++i) {
            if (i > 0) { os << ","; }
            os << v[i];
        }
        os << "]";
        return os;
    }

    /// Simple CSC-format matrix class used during simplex method
    class Matrix {
    public:
        /// Creates the 0 x 0 matrix
        Matrix() : Matrix(0,0) {}

        /// Creates a m x n all-zeros matrix 
        Matrix(size_t m, size_t n) : n_rows(m), n_cols(n) {
            for (size_t i = 0; i < n+1; ++i) {
                col_index.push_back(0);
            }
        }

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
            for (size_t c = 0; c < a.cols(); ++c) {
                for (auto it = a.begin(c); it != a.end(c); ++it) {
                    aug((*it).row, c) = (*it).value;
                }
            }
            
            for (size_t c = 0; c < b.cols(); ++c) {
                for (auto it = b.begin(c); it != b.end(c); ++it) {
                    aug((*it).row, c + a.cols()) = (*it).value;
                }
            }

            return aug;
        }

        /// Resizes the matrix. Any non-zero entries outside of the new dimensions are dropped.
        void resize(size_t r, size_t c) {
            if (c < cols() || r < rows()) { // remake the matrix
                std::vector<Number> new_v; 
                std::vector<size_t> new_cols;
                std::vector<size_t> new_rows;

                for (size_t j = 0; j < c; ++j) {
                    size_t vstart = col_index[j];
                    size_t vend = col_index[j+1];
                    new_cols.push_back(new_v.size());
                    for (size_t i = vstart; i < vend; ++i) {
                        if (row_index[i] < r) {
                            new_v.push_back(value[i]);
                            new_rows.push_back(row_index[i]);
                        }
                    }
                }
                new_cols.push_back(new_v.size());

                value = new_v;
                row_index = new_rows;
                col_index = new_cols;
            }

            while (col_index.size() < c+1) {
                col_index.push_back(value.size());
            }
        
            n_rows = r;
            n_cols = c;
        }

        /// Returns the number of rows
        size_t rows() const { return n_rows; }

        /// Returns the number of columns
        size_t cols() const { return n_cols; }

        /// Checks if the matrix is square
        bool square() const { return rows() == cols(); }

        /// Gets a single column as a vector
        Vector column(size_t c) const {
            if (c >= cols()) { throw std::out_of_range("column out-of-bounds."); }

            Vector ret(rows());
            size_t col_start = col_index[c];
            size_t col_end = col_index[c+1];
            for (size_t i = col_start; i < col_end; ++i) {
                ret[row_index[i]] = value[i];
            }
            return ret;
        }

        /// Returns a submatrix formed by a list of columns
        Matrix submatrix_cols(std::vector<size_t> subcols) const {
            Matrix ret(n_rows, subcols.size());
            ret.col_index.clear();

            for (size_t j = 0; j < subcols.size(); ++j) {
                if (subcols[j] >= n_cols) { throw std::invalid_argument("column out-of-bounds."); }

                size_t col_start = col_index[subcols[j]];
                size_t col_end = col_index[subcols[j]+1];
                ret.col_index.push_back(ret.value.size());
                for (size_t i = col_start; i < col_end; ++i) {
                    ret.value.push_back(value[i]); 
                    ret.row_index.push_back(row_index[i]);
                }
            }
            ret.col_index.push_back(ret.value.size());

            return ret;
        }

        /// Computes the row-reduced Echelon form
        void rref() {
            size_t lead = 0; 
            for (size_t r = 0; r < rows(); ++r) {
                if (lead >= cols()) { break; }

                // Find pivot
                size_t i = r;
                while (std::abs((*this)(i, lead)) < Eps) {
                    i += 1;
                    if (i == rows()) {
                        i = r;
                        lead += 1;
                        if (lead == cols()) { return; }
                    }
                }

                // Swap rows i and r
                swap_rows(i, r);

                // R_r <- R_r / pivot
                scale_row(r, 1 / (*this)(r, lead));

                // For each row i!=r, R_i <- R_i - R_r * m(i,lead)
                for (size_t i = 0; i < rows(); ++i) {
                    if (i == r) { continue; }
                    add_rows(i, r, -(*this)(i, lead));
                }

                lead += 1;
            }
        }

        Matrix inverse() const {
            if (!square()) {
                throw std::runtime_error("cannot take inverse of non-square matrix.");
            }

            Matrix aug = Matrix::augment((*this), Matrix::identity(rows()));
            aug.rref();

            // Check that matrix was invertible.
            for (size_t r = 0; r < rows(); ++r) {
                for (size_t c = 0; c < cols(); ++c) {
                    bool pass = true;
                    if (r == c && std::abs(aug(r,c) - 1) > Eps) { pass = false; }
                    if (r != c && std::abs(aug(r,c)) > Eps) { pass = false; }

                    if (!pass) {
                        throw std::runtime_error("cannot take inverse of a singular matrix.");
                    }
                }
            }

            // Get the inverse from aug
            std::vector<size_t> subcols;
            for (size_t i = 0; i < rows(); ++i) {
                subcols.push_back(i + rows());
            }
            return aug.submatrix_cols(subcols); 
        }

        /// Swaps rows r1 and r2
        void swap_rows(size_t r1, size_t r2) {
            for (size_t i = 0; i < value.size(); ++i) {
                if (row_index[i] == r1) { row_index[i] = r2; }
                else if (row_index[i] == r2) { row_index[i] = r1; }
            }
        }

        /// Scales a row by `s`
        void scale_row(size_t r, Number s) {
            for (size_t i = 0; i < value.size(); ++i) {
                if (row_index[i] == r) { value[i] *= s; } 
            }
        }

        /// Performs basic ERO R1 <- R1 + s*R2
        void add_rows(size_t r1, size_t r2, Number s) {
            /*std::vector<Number> r2_vals(cols());
            size_t c = 0;
            for (size_t i = 0; i < value.size(); ++i) {
                if (i == col_index[c+1]) { c += 1; }
                if (row_index[i] == r2) { r2_vals[c] = value[i]; }
            }

            for (size_t i = 0; i < cols(); ++i) {
                if (std::abs(r2_vals[i]) < Eps) continue;
                (*this)(r1, i) += s * r2_vals[i];
            }*/
            for (size_t j = 0; j < cols(); ++j) {
                for (size_t i = col_index[j]; i < col_index[j+1]; ++i) {
                    if (row_index[i] == r2) {
                        (*this)(r1, j) += s * value[i]; 
                    }
                }
            }
        }

        struct Entry {
            size_t row;
            Number& value;
        };
        struct ConstEntry {
            size_t row;
            const Number& value;
        };

        class iterator {
        public:
            iterator(Matrix* m, size_t idx) : m(m), idx(idx) {}
            Entry operator*() const { return {  m->row_index[idx], m->value[idx] }; }
            iterator& operator++() { ++idx; return *this; }
            bool operator!=(const iterator& other) const { return idx != other.idx; }
        private:
            Matrix* m;
            size_t idx;
        };

        class const_iterator {
        public:
            const_iterator(const Matrix* m, size_t idx) : m(m), idx(idx) {}
            ConstEntry operator*() const { return {  m->row_index[idx], m->value[idx] }; }
            const_iterator& operator++() { ++idx; return *this; }
            bool operator!=(const const_iterator& other) const { return idx != other.idx; }
        private:
            const Matrix* m;
            size_t idx;
        };

        iterator begin(size_t c) { 
            if (c >= cols()) { throw std::invalid_argument("column out-of-bounds"); }
            return iterator(this, col_index[c]); 
        }
        const_iterator begin(size_t c) const { 
            if (c >= cols()) { throw std::invalid_argument("column out-of-bounds"); }
            return const_iterator(this, col_index[c]); 
        }
        iterator end(size_t c) { 
            if (c >= cols()) { throw std::invalid_argument("column out-of-bounds"); }
            return iterator(this, col_index[c+1]); 
        }
        const_iterator end(size_t c) const { 
            if (c >= cols()) { throw std::invalid_argument("column out-of-bounds"); }
            return const_iterator(this, col_index[c+1]); 
        }

        /// Access the entry at (r,c)
        Number operator()(std::size_t r, std::size_t c) const {
            if (r >= rows() || c >= cols()) {
                throw std::invalid_argument("index out-of-bounds");
            }
            
            size_t col_start = col_index[c];
            size_t col_end = col_index[c+1];
            for (size_t i = col_start; i < col_end; ++i) {
                if (row_index[i] == r) {
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

            size_t col_start = col_index[c];
            size_t col_end = col_index[c+1];
            for (size_t i = col_start; i < col_end; ++i) {
                if (row_index[i] == r) {
                    return value[i];
                }
            }

            // Insert new entry
            size_t new_idx = col_end;
            value.insert(value.begin() + new_idx, 0); 
            row_index.insert(row_index.begin() + new_idx, r);

            for (size_t i = c+1; i < cols() + 1; ++i) {
                col_index[i] += 1;
            }
            return value[new_idx];
        }

        /// Matrix*Vector multiplication
        Vector operator*(const Vector& rhs) const {
            if (cols() != rhs.size()) {
                throw std::invalid_argument("multiply: dimension mismatch");
            }

            Vector ret(rows());
            for (size_t j = 0; j < cols(); ++j) {
                Number xj = rhs[j];
                if (std::abs(xj) < Eps) { continue; }
                for (size_t k = col_index[j]; k < col_index[j+1]; ++k) {
                    ret[row_index[k]] += value[k] * xj;
                }
            }

            return ret;
        }

        /// Compute y=x^TA
        friend Vector operator*(const Vector& lhs, const Matrix& rhs) {
            if (rhs.rows() != lhs.size()) {
                throw std::invalid_argument("multiply: dimension mismatch");
            }

            Vector ret(rhs.cols());

            for (size_t j = 0; j < rhs.cols(); ++j) {
                for (size_t k = rhs.col_index[j]; k < rhs.col_index[j+1]; ++k) {
                    ret[j] += rhs.value[k] * lhs[rhs.row_index[k]];
                }
            }

            return ret;
        }

        void operator+=(const Matrix& rhs) {
            if (cols() != rhs.cols() || rows() != rhs.rows()) {
                throw std::invalid_argument("add: dimension mismatch");
            }

            for (size_t c = 0; c < rhs.cols(); ++c) {
                for (auto it = rhs.begin(c); it != rhs.end(c); ++it) {
                    (*this)((*it).row, c) += (*it).value;
                }
            }
        }

        friend Matrix operator+(Matrix lhs, const Matrix& rhs) {
            lhs += rhs;
            return lhs;
        }

        void operator-=(const Matrix& rhs) {
            if (cols() != rhs.cols() || rows() != rhs.rows()) {
                throw std::invalid_argument("subtract: dimension mismatch");
            }

            for (size_t c = 0; c < rhs.cols(); ++c) {
                for (auto it = rhs.begin(c); it != rhs.end(c); ++it) {
                    (*this)((*it).row, c) -= (*it).value;
                }
            }
        }

        friend Matrix operator-(Matrix lhs, const Matrix& rhs) {
            lhs -= rhs;
            return lhs;
        }

        void operator*=(const Number& rhs) {
            for (size_t i = 0; i < value.size(); ++i) {
                value[i] *= rhs;
            }
        }

        friend Matrix operator*(Matrix lhs, const Number& rhs) {
            lhs *= rhs;
            return lhs;
        }
        friend Matrix operator*(Number& lhs, Matrix rhs) { return rhs * lhs; }
    private:
        std::size_t n_rows, n_cols;
        std::vector<Number> value;
        std::vector<size_t> col_index;
        std::vector<size_t> row_index;
    };

    inline Matrix operator*(const Number& lhs, const Matrix& rhs) {
        return rhs * lhs;
    }

    inline std::ostream& operator<<(std::ostream& os, const Matrix& m) {
        size_t min_w = 12;

        for (size_t r = 0; r < m.rows(); ++r) {
            for (size_t c = 0; c < m.cols(); ++c) {
                os << std::setw(min_w) << m(r, c);
            }
            os << std::endl;
        }

        return os;
    }

    enum class SolutionStatus { kOptimal, kUnbounded, kInfeasible, kFeasible, kAborted };
    enum class ProblemType { Min, Max };

    inline std::ostream& operator<<(std::ostream& os, SolutionStatus s) {
        switch (s) {
            case lp::SolutionStatus::kOptimal: return os << "Optimal";
            case lp::SolutionStatus::kUnbounded: return os << "Unbounded";
            case lp::SolutionStatus::kInfeasible: return os << "Infeasible";
            case lp::SolutionStatus::kFeasible: return os << "Feasible";
            case lp::SolutionStatus::kAborted: return os << "Aborted";
        }
        return os << "Unknown";
    }

    /// Base class solver that takes in a problem of the form:
    /// min z = c^Tx
    /// s.t. Ax == b
    ///       x >= 0
    class Solver {
    public:
        struct Solution {
            Vector x;
            Number z;
            SolutionStatus status;
        };

        virtual ~Solver() {}
        virtual Solution solve(Matrix A, Vector b, Vector c) = 0;
    };

    /// Performs revised simplex method to solve the LP.
    class DefaultSolver : public Solver {
    public:
        DefaultSolver() {} 

        Solution solve(Matrix _A, Vector _b, Vector _c) {
            A = std::move(_A);
            b = std::move(_b);
            c = _c;
            x = Vector(A.cols());

            original_num_cols = A.cols();
            status = SolutionStatus::kFeasible;
            iter_num = 0;
            start();

            while (status == SolutionStatus::kFeasible) {
                iter_num += 1;
                step();
            }

            if (cur_phase == 1) {
                Number z = c.dot(x);
                if (std::abs(z) > Eps) {
#ifdef LP_H_DEBUG
                std::cout << "Auxiliary LP solved. Problem is infeasible (z*=" << z << ")." << std::endl;
#endif
                    Vector x_soln(original_num_cols);
                    for (size_t i = 0; i < original_num_cols; ++i) { x_soln[i] = x[i]; }
                    return Solution { x_soln, z, SolutionStatus::kInfeasible };
                }
#ifdef LP_H_DEBUG
                std::cout << "Auxiliary LP solved. Problem is feasible (z*=" << z << ")." << std::endl;
#endif
                
                // TODO: Find basis without auxiliary variables.

                A.resize(A.rows(), original_num_cols);
                c = _c;
                x.resize(original_num_cols);

                // Remove auxiliary variables from bv & nbv
                for (size_t i = 0; i < bv.size(); ++i) {
                    if (bv[i] >= original_num_cols) {
                        bv.erase(bv.begin() + i);
                        i -= 1;
                    }
                }
                for (size_t i = 0; i < nbv.size(); ++i) {
                    if (nbv[i] >= original_num_cols) {
                        nbv.erase(nbv.begin() + i);
                        i -= 1;
                    }
                }

                cur_phase = 2;
                status = SolutionStatus::kFeasible;
                while (status == SolutionStatus::kFeasible) {
                    iter_num += 1;
                    step();
                }
            }

            // Create the solution
            Number z = obj_value();

#ifdef LP_H_DEBUG
            std::cout << "Final LP solved. Problem is " << status << " (z*=" << z << ")" << std::endl;
#endif

            Vector x_soln(original_num_cols);
            for (size_t i = 0; i < original_num_cols; ++i) { x_soln[i] = x[i]; }
            return Solution { x_soln, z, status };
        }
    protected:
        Number obj_value() const {
            return c.dot(x);
        }

        virtual void start() {
#ifdef LP_H_DEBUG
            std::cout << "c = " << c << std::endl;
            std::cout << "A = " << std::endl << A;
            std::cout << "b = " << b << std::endl;
#endif
            std::vector<size_t> identity_cols(A.rows(), false);
            std::vector<bool> is_identity_cols(A.rows(), false);

            // Quickly check if the identity matrix is a submatrix, if so, set those columns as bv
            /*if (A.rows() > 1) {
                for (size_t i = 0; i < A.cols(); ++i) {
                    Number nonzero_entry = 0;
                    size_t num_zeros = A.rows();
                    size_t nonzero_row = 0;
                    for (auto it = A.begin(i); it != A.end(i); ++it) {
                        if (std::abs((*it).value) > Eps) { 
                            num_zeros -= 1;
                            nonzero_entry = (*it).value; 
                            nonzero_row = (*it).row;
                        }
                    }

                    if (num_zeros == A.rows()-1 && std::abs(nonzero_entry-1) < Eps && bv.size() < A.rows()) {
                        is_identity_cols[nonzero_row] = true;
                        identity_cols[nonzero_row] = i;
                    }
                }
            }*/

            // Add artificial variables
            if (true) { //(std::find(is_identity_cols.begin(), is_identity_cols.end(), false) != is_identity_cols.end()) {
                cur_phase = 1;
                
#ifdef LP_H_DEBUG
                std::cout << "No identity matrix found, adding artificial variables." << std::endl;
#endif
                A.resize(A.rows(), A.cols() + A.rows());
                x.resize(c.size() + A.rows());
                c = Vector(c.size() + A.rows());
                for (size_t i = 0; i < A.rows(); ++i) {
                    bv.push_back(original_num_cols + i);
                    A(i, original_num_cols + i) = 1;
                    c[original_num_cols + i] = 1;
                }

                for (size_t i = 0; i < original_num_cols; ++i) { nbv.push_back(i); }
#ifdef LP_H_DEBUG
                std::cout << "New problem:" << std::endl;
                std::cout << "c = " << c << std::endl;
                std::cout << "A = " << std::endl << A;
                std::cout << "b = " << b << std::endl;
#endif
            }
            else {
                for (size_t i = 0; i < A.rows(); ++i) {
                    bv.push_back(identity_cols[i]);
                }
                for (size_t i = 0; i < A.rows(); ++i) {
                    if (std::find(bv.begin(), bv.end(), i) == bv.end()) {
                        nbv.push_back(i);
                    }
                }
                cur_phase = 2;
            }

            // Compute basis matrix and invert
            Binv = A.submatrix_cols(bv);

            // Determine initial bfs
            for (size_t i = 0; i < b.size(); ++i) {
                x[bv[i]] = b[i];
            }
        }

        virtual void step() {
#ifdef LP_H_DEBUG
            // Print the current status
            std::cout << std::endl << "Iteration #" << iter_num << ":" << std::endl;
            std::cout << "----------------------" << std::endl;
            std::cout << "bv = [";
            for (size_t i = 0; i < bv.size(); ++i) { 
                if (i > 0) { std::cout << ","; }
                std::cout << bv[i];
            }
            std::cout << "]" << std::endl << "nbv = [";
            for (size_t i = 0; i < nbv.size(); ++i) { 
                if (i > 0) { std::cout << ","; }
                std::cout << nbv[i];
            }
            std::cout << "]" << std::endl << "x = " << x;
            std::cout << std::endl << "Binv = " << std::endl << Binv;
#endif

            // Find the new objective vector
            Vector cB = c.subvector(bv);
            Vector cP = Vector(nbv.size());
            Vector y = cB * Binv;
#ifdef LP_H_DEBUG
            std::cout << "y=" << y << std::endl;
#endif

            for (size_t j = 0; j < nbv.size(); ++j) {
                cP[j] = c[nbv[j]] - y.dot(A.column(nbv[j]));
            }

#ifdef LP_H_DEBUG
            std::cout << "c' = " << cP << std::endl;
#endif

            // Check if optimal
            size_t entering_var = std::numeric_limits<size_t>::max();
            for (size_t i = 0; i < cP.size(); ++i) {
                if (cP[i] < -Eps) { 
                    if (nbv[i] < entering_var) {
                        entering_var = nbv[i];
                    }
                }
            }

            if (entering_var == std::numeric_limits<size_t>::max()) {
                status = SolutionStatus::kOptimal;
                return;
            }

#ifdef LP_H_DEBUG
            std::cout << "entering_var = " << entering_var << std::endl;
#endif

            // Compute feasible direction vector and check if unbounded.
            Vector d(x.size());
            Vector dB = -1.0 * Binv * A.column(entering_var);
            bool unbounded = true;
            for (size_t i = 0; i < dB.size(); ++i) {
                if (dB[i] < -Eps) { unbounded = false; }
                d[bv[i]] = dB[i];
            }

            if (unbounded) {
                status = SolutionStatus::kUnbounded;
                return;
            }
            d[entering_var] = 1;

#ifdef LP_H_DEBUG
            std::cout << "d = " << d << std::endl;
#endif

            // Perform minimum-ratio test to find leaving variable
            size_t leaving_var = std::numeric_limits<size_t>::max();
            Number theta = Infinity;
            for (size_t i = 0; i < dB.size(); ++i) {
                if (dB[i] >= -Eps) { continue; }
                Number r = -x[bv[i]] / dB[i];

#ifdef LP_H_DEBUG
                std::cout << "i=" << bv[i] <<": " << -x[bv[i]] << "/" << dB[i] << "=" << r << std::endl;
#endif

                if (r < theta - Eps || (std::abs(r - theta) < Eps && bv[i] < bv[leaving_var])) {
                    theta = r;
                    leaving_var = i;
                }
            }

            if (leaving_var == std::numeric_limits<size_t>::max()) {
                status = SolutionStatus::kUnbounded;
                return;
            }
#ifdef LP_H_DEBUG
            std::cout << "leaving_var = " << bv[leaving_var] << std::endl;
            std::cout << "theta = " << theta << std::endl;
#endif

            if (theta < -Eps) { throw std::logic_error("negative theta detected."); }

            // Compute the new bfs
            x += theta * d;
            x[bv[leaving_var]] = 0;

            // Update Binv using EROs
            Binv.scale_row(leaving_var, -1/dB[leaving_var]);
            for (size_t r = 0; r < Binv.rows(); ++r) {
                if (r == leaving_var) { continue; } 
                Binv.add_rows(r, leaving_var, dB[r]);
            }

            // Finally update bv and nbv
            size_t leaving_var_actual = bv[leaving_var];
            bv[leaving_var] = entering_var;
            for (size_t i = 0; i < nbv.size(); ++i) {
                if (nbv[i] == entering_var) { nbv[i] = leaving_var_actual; break; }
            }
        }

    protected:
        Matrix A;
        Vector b;
        Vector c;
        Vector x;
        size_t iter_num;
        std::vector<size_t> bv;
        std::vector<size_t> nbv;
        Matrix Binv;
        SolutionStatus status;
        size_t cur_phase;
    private:
        size_t original_num_cols;
    };

    struct Variable {
        size_t id;
        Number min;
        Number max;
        std::string name;    
    };
    
    inline std::ostream& operator<<(std::ostream& os, const Variable& v) {
        os << v.name;
        return os;
    }

    struct Expression {
        std::vector<std::pair<Number, Variable*>> terms;

        Expression() {}

        Expression& operator+=(Expression rhs) {
            for (auto& rhs_term : rhs.terms) {
                auto it = std::find_if(terms.begin(), terms.end(), [&](std::pair<Number, Variable*> term) {
                    return term.second->id == rhs_term.second->id;
                });
                if (it != terms.end()) {
                    it->first += rhs_term.first;
                }
                else {
                    terms.push_back(rhs_term);         
                }
            }

            return *this;
        }
        
        Expression& operator-=(Expression rhs) {
            for (auto& rhs_term : rhs.terms) {
                auto it = std::find_if(terms.begin(), terms.end(), [&](std::pair<Number, Variable*> term) {
                    return term.second->id == rhs_term.second->id;
                });
                if (it != terms.end()) {
                    it->first -= rhs_term.first;
                }
                else {
                    terms.push_back(std::make_pair(-rhs_term.first, rhs_term.second));         
                }
            }

            return *this;
        }

        Expression& operator*=(Number rhs) {
            if (rhs == 0) { 
                terms.clear();
            }
            
            for (auto& term: terms) {
                term.first *= rhs;
            }

            return *this;
        }

        Expression& operator+=(Variable& rhs) {
            auto it = std::find_if(terms.begin(), terms.end(), [&](std::pair<Number, Variable*> term) {
                return term.second->id == rhs.id;
            });
            if (it != terms.end()) {
                it->first += 1;
            }
            else {
                terms.push_back(std::make_pair(1, &rhs));         
            }

            return *this;
        }

        Expression& operator-=(Variable& rhs) {
            auto it = std::find_if(terms.begin(), terms.end(), [&](std::pair<Number, Variable*> term) {
                return term.second->id == rhs.id;
            });
            if (it != terms.end()) {
                it->first -= 1;
            }
            else {
                terms.push_back(std::make_pair(-1, &rhs));         
            }

            return *this;
        }
        
        Expression operator-() { return (*this) * Number(-1); }
        friend Expression operator+(Variable& lhs, Expression rhs) { return rhs + lhs; }
        friend Expression operator+(Expression lhs, Variable& rhs) {
            lhs += rhs;
            return lhs;
        }
        friend Expression operator-(Variable& lhs, Expression rhs) { return rhs - lhs; }
        friend Expression operator-(Expression lhs, Variable& rhs) {
            lhs -= rhs;
            return lhs;
        }
        friend Expression operator+(Expression lhs, Expression rhs) {
            lhs += rhs;
            return lhs;
        }
        friend Expression operator-(Expression lhs, Expression rhs) {
            lhs -= rhs;
            return lhs;
        }
        friend Expression operator*(Number lhs, Expression rhs) {
            rhs *= lhs;
            return rhs;
        }
        friend Expression operator*(Expression lhs, Number rhs) { return rhs * lhs; }
    };

    inline Expression operator*(Variable& lhs, Number rhs) {
        Expression e; 
        e.terms.push_back(std::make_pair(rhs, &lhs));
        return e;
    }
    inline Expression operator*(Number lhs, Variable& rhs) { return rhs * lhs;} 
    inline Expression operator-(Variable& v) { return v * Number(-1); }
    inline Expression operator+(Variable& lhs, Variable& rhs) { 
        Expression e;
        e += lhs;
        e += rhs;
        return e;
    } 
    inline Expression operator-(Variable& lhs, Variable& rhs) { 
        Expression e;
        e += lhs;
        e -= rhs;
        return e;
    } 

    inline std::ostream& operator<<(std::ostream& os, const Expression& e) {
        bool first_term = true;

        if (e.terms.size() == 0) { os << "0"; }

        for (size_t i = 0; i < e.terms.size(); ++i) {
            Number c = e.terms[i].first;
            if (c == 0) { continue; }

            if (!first_term) {
                if (c > 0) { os << " + "; }        
                if (c < 0) { os << " - "; }        
            }
            else if (c < 0) { os << "-"; }
            if (std::abs(c) != 1) { os << std::abs(c) << "*"; }
            os << e.terms[i].second->name;
            first_term = false;
        }
        return os;
    }

    enum class ConstraintType { Eq, LEq, GEq };

    struct Constraint {
        Expression lhs;
        Number rhs;
        ConstraintType type;
    };

    inline Constraint operator==(Expression lhs, Number rhs) {
        Constraint c { std::move(lhs), rhs, ConstraintType::Eq };
        return c;
    }
    inline Constraint operator<=(Expression lhs, Number rhs) {
        Constraint c { std::move(lhs), rhs, ConstraintType::LEq };
        return c;
    }
    inline Constraint operator>=(Expression lhs, Number rhs) {
        Constraint c { std::move(lhs), rhs, ConstraintType::GEq };
        return c;
    }
    inline Constraint operator==(Variable& lhs, Number rhs) {
        Constraint c { 1*lhs, rhs, ConstraintType::Eq };
        return c;
    }
    inline Constraint operator<=(Variable& lhs, Number rhs) {
        Constraint c { 1*lhs, rhs, ConstraintType::LEq };
        return c;
    }
    inline Constraint operator>=(Variable& lhs, Number rhs) {
        Constraint c { 1*lhs, rhs, ConstraintType::GEq };
        return c;
    }

    inline std::ostream& operator<<(std::ostream& os, const Constraint& c) {
        os << c.lhs << " ";
        switch (c.type) {
            case ConstraintType::Eq: os << "="; break;
            case ConstraintType::LEq: os << "<="; break;
            case ConstraintType::GEq: os << ">="; break;
        }
        os << " " << c.rhs;
        return os;
    }


    struct Solution {
        Number objective;
        SolutionStatus status;
        std::vector<Number> values;
    
        Number operator[](const Variable& v) {
            if (v.id >= values.size()) { 
                throw std::runtime_error("Variable not part of problem.");
            }
            return values[v.id];
        }
    };

    namespace transformations {
        class Base {
        public:
            Base(bool additional_var) : additional_var(additional_var) {}
            virtual ~Base() {}
            virtual void apply_matrix(Matrix&) {}
            virtual void apply_rhs(Vector&) {}
            virtual void apply_obj(Vector&) {}
            virtual void apply_soln(Solver::Solution&) {}
            virtual bool adds_var() { return additional_var; }
        protected:
            bool additional_var;
        };

        class NegateVar : public Base {
        public:
            NegateVar(size_t id) : Base(false), id(id) {}
            void apply_matrix(Matrix& A) override {
                for (size_t r = 0; r < A.rows(); ++r) {
                    A(r, id) *= -1;
                }
            }
    
            void apply_obj(Vector& c) override {
                c[id] *= -1;             
            }

            void apply_soln(Solver::Solution& s) override {
                s.x[id] *= -1;             
            }
        private:
            size_t id;
        };

        class NegateRow : public Base {
        public:
            NegateRow(size_t row) : Base(false), row(row) {}
            void apply_matrix(Matrix& A) override {
                A.scale_row(row, -1);
            }
            
            void apply_rhs(Vector& b) override {
                b[row] *= -1;  
            }
        private:
            size_t row;
        };

        class VariableBound : public Base {
        public:
            VariableBound(size_t constraint_id, size_t var_id, size_t slack_var, ConstraintType type, Number bound) : Base(true), 
                constraint_id(constraint_id), var_id(var_id), slack_var(slack_var), type(type), bound(bound) {}

            void apply_matrix(Matrix& A) override {
                A(constraint_id, var_id) = 1;
                A(constraint_id, slack_var) = type == ConstraintType::LEq ? 1 : -1;
            };

            void apply_rhs(Vector& b) override {
                b[constraint_id] = bound;
            };
        private:
            size_t constraint_id;
            size_t var_id;
            size_t slack_var;
            ConstraintType type;
            Number bound;
        };

        class SlackVar : public Base {
        public:
            SlackVar(size_t id, size_t row, ConstraintType type) : Base(true), id(id), row(row), type(type) {}
            void apply_matrix(Matrix& A) override {
                A(row, id) = type == ConstraintType::LEq ? 1 : -1;
            }
        private:
            size_t id;
            size_t row;
            ConstraintType type;
        };

        class URSVar : public Base {
        public:
            URSVar(size_t vplus, size_t vneg) : Base(true), vplus(vplus), vneg(vneg) {}
            void apply_matrix(Matrix& A) override {
                for (size_t r = 0; r < A.rows(); ++r) {
                    Number v = A(r, vplus);
                    if (std::abs(v) > Eps) {
                        A(r, vneg) = -v; 
                    }
                }
            }

            void apply_obj(Vector& c) override {
                Number v = c[vplus];
                if (std::abs(v) > Eps) {
                    c[vneg] = -v;
                }
            }

            void apply_soln(Solver::Solution& s) override {
                s.x[vplus] -= s.x[vneg];             
            }
        private:
            size_t vplus;
            size_t vneg;
        };
    }

    class Problem {
    public:
        Problem(ProblemType type) : type(type) {}

        static Problem minimize() { return Problem(ProblemType::Min); }
        static Problem maximize() { return Problem(ProblemType::Max); }

        Variable add_var(Number min, Number max) {
            size_t id = variables.size();
            return add_var(min, max, "x" + std::to_string(id));
        }

        Variable add_var(Number min, Number max, std::string name) {
            size_t id = variables.size();
            variables.push_back(Variable { id, min, max, name });
            return variables[id];
        }

        void objective(Expression e) {
            objective_ = std::move(e);
        }

        void constraint(Constraint c) {
            constraints.push_back(std::move(c));
        }

        Solution solve() {
            DefaultSolver * solver = new DefaultSolver();
            Solution soln = solve(solver);
            delete solver;
            return soln;
        }

        Solution solve(Solver * solver) {
            // ------------------------
            // Convert to standard form
            // ------------------------
            std::vector<transformations::Base*> transformations;
            size_t num_extra_vars = 0;
            size_t num_extra_constraints = 0;
            
            // 1. Convert non-equality constraints to equality with slack variables
            for (size_t r = 0; r < constraints.size(); ++r) {
                if (constraints[r].type != ConstraintType::Eq) {
#ifdef LP_H_DEBUG
                    std::cout << "Adding slack variable x" << variables.size() + num_extra_vars << std::endl;
#endif
                    transformations.push_back(new transformations::SlackVar(variables.size() + num_extra_vars, r, constraints[r].type));
                    num_extra_vars += 1;
                }
            }

            // 2. Variable bounds
            for (size_t c = 0; c < variables.size(); ++c) {
                if (std::abs(variables[c].max) < Eps) { // Variable is negative
#ifdef LP_H_DEBUG
                    std::cout << "Variable x" << c << " is negative, negating." << std::endl;
#endif
                    transformations.push_back(new transformations::NegateVar(c));

                    if (variables[c].min > -Infinity) {
#ifdef LP_H_DEBUG
                    std::cout << "Variable x" << c << " is bounded below, adding slack var x" << variables.size() + num_extra_vars << " and constraint c" << constraints.size() + num_extra_constraints << std::endl;
#endif
                        transformations.push_back(new transformations::VariableBound(constraints.size() + num_extra_constraints,
                                    c, variables.size() + num_extra_vars, ConstraintType::LEq, -variables[c].min));            
                        transformations.push_back(new transformations::NegateRow(constraints.size() + num_extra_constraints));
                        num_extra_constraints += 1;
                        num_extra_vars += 1;
                    }
                }
                else if (variables[c].min < -Eps) { // Variable is URS
#ifdef LP_H_DEBUG
                    std::cout << "Variable x" << c << " is urs, adding negative x" << variables.size() + num_extra_vars << std::endl;
#endif
                    transformations.push_back(new transformations::URSVar(c, variables.size() + num_extra_vars));
                    size_t cneg = variables.size() + num_extra_vars;
                    num_extra_vars += 1;

                    if (variables[c].min > -Infinity) {
#ifdef LP_H_DEBUG
                    std::cout << "Variable x" << c << " is bounded below, adding slack var x" << variables.size() + num_extra_vars << " and constraint c" << constraints.size() + num_extra_constraints << std::endl;
#endif
                        transformations.push_back(new transformations::VariableBound(constraints.size() + num_extra_constraints,
                                    c, variables.size() + num_extra_vars, ConstraintType::GEq, variables[c].min));            
#ifdef LP_H_DEBUG
                    std::cout << "Variable x" << cneg << " is bounded below, adding slack var x" << variables.size() + num_extra_vars + 1 << " and constraint c" << constraints.size() + num_extra_constraints + 1 << std::endl;
#endif
                        transformations.push_back(new transformations::VariableBound(constraints.size() + num_extra_constraints + 1,
                                    cneg, variables.size() + num_extra_vars + 1, ConstraintType::GEq, variables[c].min));            
                        transformations.push_back(new transformations::NegateRow(constraints.size() + num_extra_constraints));
                        transformations.push_back(new transformations::NegateRow(constraints.size() + num_extra_constraints + 1));
                        num_extra_constraints += 2;
                        num_extra_vars += 2;
                    }

                    if (variables[c].max < Infinity) {
#ifdef LP_H_DEBUG
                    std::cout << "Variable x" << c << " is bounded above, adding slack var x" << variables.size() + num_extra_vars << " and constraint c" << constraints.size() + num_extra_constraints << std::endl;
#endif
                        transformations.push_back(new transformations::VariableBound(constraints.size() + num_extra_constraints,
                                    c, variables.size() + num_extra_vars, ConstraintType::GEq, variables[c].min));            
#ifdef LP_H_DEBUG
                    std::cout << "Variable x" << cneg << " is bounded above, adding slack var x" << variables.size() + num_extra_vars + 1 << " and constraint c" << constraints.size() + num_extra_constraints + 1 << std::endl;
#endif
                        transformations.push_back(new transformations::VariableBound(constraints.size() + num_extra_constraints,
                                    c, variables.size() + num_extra_vars, ConstraintType::LEq, variables[c].max));            
                        transformations.push_back(new transformations::VariableBound(constraints.size() + num_extra_constraints + 1,
                                    cneg, variables.size() + num_extra_vars + 1, ConstraintType::LEq, variables[c].max));            
                        num_extra_constraints += 2;
                        num_extra_vars += 2;
                    }
                }
                else { // Variable is positive
                    if (variables[c].max < Infinity) {
#ifdef LP_H_DEBUG
                    std::cout << "Variable x" << c << " is bounded above, adding slack var x" << variables.size() + num_extra_vars << " and constraint c" << constraints.size() + num_extra_constraints << std::endl;
#endif
                        transformations.push_back(new transformations::VariableBound(constraints.size() + num_extra_constraints,
                                    c, variables.size() + num_extra_vars, ConstraintType::LEq, variables[c].max));            
                        num_extra_constraints += 1;
                        num_extra_vars += 1;
                    }
                    if (variables[c].min > Eps) {
#ifdef LP_H_DEBUG
                    std::cout << "Variable x" << c << " is bounded below, adding slack var x" << variables.size() + num_extra_vars << " and constraint c" << constraints.size() + num_extra_constraints << std::endl;
#endif
                        transformations.push_back(new transformations::VariableBound(constraints.size() + num_extra_constraints,
                                    c, variables.size() + num_extra_vars, ConstraintType::GEq, variables[c].min));            
                        num_extra_constraints += 1;
                        num_extra_vars += 1;
                    }
                }
            }

            // 3. Ensure b is positive
            for (size_t i = 0; i < constraints.size(); ++i) {
                if (constraints[i].rhs < -Eps) {
                    transformations.push_back(new transformations::NegateRow(i));
                }
            }
            
            // ---------------------
            // Build the matrix form
            // ---------------------
            Matrix A(constraints.size() + num_extra_constraints, variables.size() + num_extra_vars);
            Vector b(constraints.size() + num_extra_constraints);
            for (size_t i = 0; i < constraints.size(); ++i) {
                for (const auto& v : constraints[i].lhs.terms) {
                    if (v.first == 0) { continue; }
                    A(i, v.second->id) = v.first;
                }

                b[i] = constraints[i].rhs;
            }

            
            Vector c(variables.size() + num_extra_vars);
            for (const auto& v: objective_.terms) {
                if (v.first == 0) { continue; }
                c[v.second->id] = (type == ProblemType::Max) ? -v.first :  v.first;
            }
#ifdef LP_H_DEBUG
            std::cout << "Problem before transformations:" << std::endl;
            std::cout << "c = " << c << std::endl ;
            std::cout << "A = " << std::endl << A;
            std::cout << "b = " << b << std::endl ;
#endif

            // Apply transformations
            for (auto& t : transformations) {
                t->apply_matrix(A);
                t->apply_rhs(b);
                t->apply_obj(c);
            }

            // -----
            // Solve
            // -----
            Solver::Solution s = solver->solve(std::move(A), std::move(b), std::move(c));

            for (auto& t : transformations) {
                t->apply_soln(s);
                delete t;
            }

            Solution solution;
            solution.status = s.status;
            for (size_t i = 0; i < variables.size(); ++i) {
                solution.values.push_back(s.x[i]);
            }
            solution.objective = (type == ProblemType::Max) ? -s.z : s.z;
            return solution;
        }

        friend std::ostream& operator<<(std::ostream& os, const Problem& p) {
            switch (p.type) {
                case ProblemType::Min: os << "min"; break;
                case ProblemType::Max: os << "max"; break;
            }
            os << "  " << p.objective_;
            os << std::endl << "s.t. ";
            for (size_t i = 0; i < p.constraints.size(); ++i) {
                if (i > 0) { os << "     "; }
                os << p.constraints[i] << std::endl;
            }
            for (size_t i = 0; i < p.variables.size(); ++i) {
                if (p.variables[i].min == -Infinity && p.variables[i].max == Infinity) { continue; }
                
                os << "     "; 
                if (p.variables[i].min == -Infinity) { os << p.variables[i] << " <= " << p.variables[i].max; }
                else if (p.variables[i].max == Infinity) { os << p.variables[i] << " >= " << p.variables[i].min; }
                else { os << p.variables[i].min << " <= " << p.variables[i] << " <= " << p.variables[i].max; }
                os << std::endl;
            }

            return os;
        }
    private:
        ProblemType type;
        Expression objective_;
        std::vector<Variable> variables;
        std::vector<Constraint> constraints;
    };
}

