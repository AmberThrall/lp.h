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

    /// Simple COO-format matrix class used during simplex method
    class Matrix {
    public:
        /// Creates the 0 x 0 matrix
        Matrix() : n_rows(0), n_cols(0) {}

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

        /// Resizes the matrix. Any non-zero entries outside of the new dimensions are dropped.
        void resize(size_t r, size_t c) {
            for (size_t i = 0; i < value.size(); ++i) {
                if (row[i] >= r || col[i] >= c) {
                    value.erase(value.begin() + i);
                    col.erase(col.begin() + i);
                    row.erase(row.begin() + i);
                    i -= 1;
                }
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
            Entry operator*() const { return {  m->row[idx], m->col[idx], m->value[idx] }; }
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
            std::vector<Number> x;
            Number z;
            SolutionStatus status;
        };

        virtual ~Solver() {}
        virtual Solution solve(Matrix A, Matrix b, Matrix c) = 0;
    };

    /// Performs revised simplex method to solve the LP.
    class DefaultSolver : public Solver {
    public:
        DefaultSolver() {} 

        Solution solve(Matrix _A, Matrix _b, Matrix _c) {
            A = std::move(_A);
            b = std::move(_b);
            c = _c;
            x = Matrix(A.cols(), 1);
            original_num_cols = A.cols();

            status = SolutionStatus::kFeasible;
            iter_num = 0;
            start();
            cP = c;

            while (status == SolutionStatus::kFeasible) {
                iter_num += 1;
                step();
            }

            if (A.cols() > original_num_cols) {
                Number z = 0;
                for (size_t i = 0; i < cP.cols(); ++i) {
                    z += cP(0, i) * x(i, 0);
                }
                if (z > Eps) {
#ifdef LP_H_DEBUG
                std::cout << "Auxiliary LP solved. Problem is infeasible (z*=" << z << ")." << std::endl;
#endif
                    std::vector<Number> x_soln(original_num_cols);
                    for (const auto& e: x) { x_soln[e.row] = e.value; }
                    return Solution { x_soln, z, SolutionStatus::kInfeasible };
                }
#ifdef LP_H_DEBUG
                std::cout << "Auxiliary LP solved. Problem is feasible (z*=" << z << ")." << std::endl;
#endif

                A.resize(A.rows(), original_num_cols);
                c = _c;
                cP = c;
                x.resize(original_num_cols, 1);

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

                status = SolutionStatus::kFeasible;
                while (status == SolutionStatus::kFeasible) {
                    iter_num += 1;
                    step();
                }
            }

            // Check feasibility
            Number z = 0;
            for (size_t i = 0; i < cP.cols(); ++i) {
                z += cP(0, i) * x(i, 0);
            }
            if (z > Eps) { status = SolutionStatus::kInfeasible; }

            // Create the solution
            z = obj_value();
            std::vector<Number> x_soln(original_num_cols);
            for (const auto& e: x) { x_soln[e.row] = e.value; }

            return Solution { x_soln, z, status };
        }
    protected:
        Number obj_value() const {
            Number obj = 0;
            for (size_t i = 0; i < c.cols(); ++i) {
                obj += c(0, i) * x(i, 0);
            }
            return obj;
        }

        virtual void start() {
#ifdef LP_H_DEBUG
            std::cout << "c = " << std::endl << c;
            std::cout << "A = " << std::endl << A;
            std::cout << "b = " << std::endl << b;
#endif

            // TODO: Fix this identify code
            // Quickly check if the identity matrix is a submatrix, if so, set those columns as bv
            for (size_t i = 0; i < A.cols(); ++i) {
                double nonzero_entry = 0;
                size_t num_zeros = 0;
                for (size_t j = 0; j < A.rows(); ++j) {
                    if (std::abs(A(j,i)) < Eps) { 
                        num_zeros += 1;
                    }
                    else { nonzero_entry = A(j,i); }
                }

                if (num_zeros == A.rows()-1 && std::abs(nonzero_entry-1) < Eps && bv.size() < A.rows()) {
                    bv.push_back(i);
                }
                else {
                    nbv.push_back(i);
                }
            }

            // Add artificial variables
            if (bv.size() != A.rows()) {
                bv.clear();
                nbv.clear();
                
#ifdef LP_H_DEBUG
                std::cout << "No identity matrix found, adding artificial variables." << std::endl;
#endif
                A.resize(A.rows(), A.cols() + A.rows());
                c = Matrix(1, c.cols() + A.rows());
                x.resize(c.cols() + A.rows(), 1);
                for (size_t i = 0; i < A.rows(); ++i) {
                    bv.push_back(original_num_cols + i);
                    A(i, original_num_cols + i) = 1;
                    c(0, original_num_cols + i) = 1;
                }

                for (size_t i = 0; i < original_num_cols; ++i) { nbv.push_back(i); }
#ifdef LP_H_DEBUG
                std::cout << "New problem:" << std::endl;
                std::cout << "c = " << std::endl << c;
                std::cout << "A = " << std::endl << A;
                std::cout << "b = " << std::endl << b;
#endif
            }

            // Compute basis matrix and invert
            Binv = A.submatrix_cols(bv);

            // Determine initial bfs
            for (const auto & e : b) {
                x(bv[e.row],0) = e.value;
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
            std::cout << "]" << std::endl << "x = <";
            for (size_t i = 0; i < x.rows(); ++i) { 
                if (i > 0) { std::cout << ","; }
                std::cout << x(i,0);
            }
            std::cout << ">" << std::endl << "c = <";
            for (size_t i = 0; i < cP.cols(); ++i) { 
                if (i > 0) { std::cout << ","; }
                std::cout << cP(0,i);
            }
            std::cout << ">" << std::endl << "Binv = " << std::endl << Binv;
#endif

            // Find the new objective vector
            Matrix cB = cP.submatrix_cols(bv);

            cP -= cB * Binv * A;
            Matrix cN = cP.submatrix_cols(nbv);

#ifdef LP_H_DEBUG
            std::cout << std::endl << "c' = " << cP;
#endif

            // Check if optimal
            size_t entering_var = 0;
            Number most_negative = 0;
            for (const auto& e : cN) { 
                if (e.value < -Eps) { 
                    if (e.value < most_negative) {
                        entering_var = nbv[e.col];
                        most_negative = e.value;
                    }
                }
            }

            if (most_negative == 0) {
                status = SolutionStatus::kOptimal;
                return;
            }

#ifdef LP_H_DEBUG
            std::cout << "entering_var = " << entering_var << std::endl;
#endif

            // Compute feasible direction vector and check if unbounded.
            Matrix d(x.rows(), 1);
            Matrix dB = Binv * A.submatrix_cols({ entering_var }) * -1.0;
            bool unbounded = true;
            for (const auto& e : dB) {
                if (e.value < -Eps) { unbounded = false; }
                d(bv[e.row], 0) = e.value;
            }

            if (unbounded) {
                status = SolutionStatus::kUnbounded;
                return;
            }
            d(entering_var, 0) = 1;

#ifdef LP_H_DEBUG
            std::cout << "d = " << std::endl << d;
#endif

            // Perform minimum-ratio test to find leaving variable
            size_t leaving_var = 0;
            Number theta = 0;
            for (const auto& e : dB) {
                if (e.value > -Eps) { continue; }
                Number r = -x(bv[e.row], 0) / e.value;

#ifdef LP_H_DEBUG
                std::cout << "i=" << bv[e.row] <<": " << -x(bv[e.row],0) << "/" << e.value << "=" << r << std::endl;
#endif

                if (r < theta || theta == 0) {
                    theta = r;
                    leaving_var = e.row;
                }
            }

#ifdef LP_H_DEBUG
            std::cout << "leaving_var = " << bv[leaving_var] << std::endl;
            std::cout << "theta = " << theta << std::endl;
#endif

            // Compute the new bfs
            x += theta * d;

            // Update Binv using EROs
            Matrix aug = Matrix::augment(Binv, dB * -1);
            
            aug.scale_row(leaving_var, -1/dB(leaving_var, 0));
            for (size_t r = 0; r < aug.rows(); ++r) {
                if (r == leaving_var) { continue; } 
                aug.add_rows(r, leaving_var, dB(r, 0));
            }

            std::vector<size_t> subcols;
            for (size_t i = 0; i < aug.cols() - 1; ++i) { subcols.push_back(i); }
            Binv = aug.submatrix_cols(subcols);

            // Finally update bv and nbv
            size_t leaving_var_actual = bv[leaving_var];
            bv[leaving_var] = entering_var;
            for (size_t i = 0; i < nbv.size(); ++i) {
                if (nbv[i] == entering_var) { nbv[i] = leaving_var_actual; break; }
            }
        }

    protected:
        Matrix A;
        Matrix b;
        Matrix c;
        Matrix cP;
        Matrix x;
        size_t iter_num;
        std::vector<size_t> bv;
        std::vector<size_t> nbv;
        Matrix Binv;
        SolutionStatus status;
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
        friend Expression operator+(Variable& lhs, Expression& rhs) { return rhs + lhs; }
        friend Expression operator+(Expression lhs, Variable& rhs) {
            lhs += rhs;
            return lhs;
        }
        friend Expression operator-(Variable& lhs, Expression& rhs) { return rhs - lhs; }
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
            virtual void apply_obj(Matrix&) {}
            virtual void apply_soln(Solver::Solution&) {}
            virtual bool adds_var() { return additional_var; }
            virtual size_t var_id() { return 0; }
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
    
            void apply_obj(Matrix& c) override {
                c(0, id) *= -1;             
            }

            void apply_soln(Solver::Solution& s) override {
                s.x[id] *= -1;             
            }
        private:
            size_t id;
        };

        class SlackVar : public Base {
        public:
            SlackVar(size_t id, size_t row, ConstraintType type) : Base(true), id(id), row(row), type(type) {}
            void apply_matrix(Matrix& A) override {
                A(row, id) = type == ConstraintType::LEq ? 1 : -1;
            }

            size_t var_id() override { return id; }
        private:
            size_t id;
            size_t row;
            ConstraintType type;
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
            
            // 1. Convert non-equality constraints to equality with slack variables
            for (size_t r = 0; r < constraints.size(); ++r) {
                if (constraints[r].type != ConstraintType::Eq) {
                    transformations.push_back(new transformations::SlackVar(variables.size() + num_extra_vars, r, constraints[r].type));
                    num_extra_vars += 1;
                }
            }

            // 2. Negate negative variables.
            for (size_t c = 0; c < variables.size(); ++c) {
                if (variables[c].max == 0.0) { 
                    transformations.push_back(new transformations::NegateVar(c));
                }
            }
            
            // Build the matrix form
            Matrix A(constraints.size(), variables.size() + num_extra_vars);
            Matrix b(constraints.size(), 1);
            for (size_t i = 0; i < constraints.size(); ++i) {
                for (const auto& v : constraints[i].lhs.terms) {
                    if (v.first == 0) { continue; }
                    A(i, v.second->id) = v.first;
                }

                b(i,0) = constraints[i].rhs;
            }

            
            Matrix c(1, variables.size() + num_extra_vars);
            for (const auto& v: objective_.terms) {
                if (v.first == 0) { continue; }
                c(0, v.second->id) = (type == ProblemType::Max) ? -v.first :  v.first;
            }

            for (auto& t : transformations) {
                t->apply_matrix(A);
                t->apply_obj(c);
            }

            // Solve
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

