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
#include <limits>

namespace lp {
#ifdef LP_H_USE_FLOAT
    using Number = float;
#else
    using Number = double;
#endif
    constexpr Number Inf = std::numeric_limits<Number>::infinity();

    /// Simple COO-format matrix class used during simplex method
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

    inline std::ostream& operator<<(std::ostream& os, const Matrix& m) {
        size_t min_w = 8;

        for (size_t r = 0; r < m.rows(); ++r) {
            for (size_t c = 0; c < m.cols(); ++c) {
                os << std::setw(min_w) << m(r, c);
            }
            os << std::endl;
        }

        return os;
    }

    enum class SolutionStatus { kOptimal, kUnbounded, kInfeasible, kFeasible, kAborted };

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

            for (size_t i = 0; i < state.nbv.size(); ++i) {
                if (state.nbv[i] == enter_var) { state.nbv[i] = leaving_var; break; }
            }

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

    struct Expression;

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

    enum class ProblemType { Min, Max };

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
            // Convert to standard form
            
            // Build the matrix form
            Matrix A(constraints.size(), variables.size());
            Matrix b(constraints.size(), 1);
            for (size_t i = 0; i < constraints.size(); ++i) {
                for (const auto& v : constraints[i].lhs.terms) {
                    if (v.first == 0) { continue; }
                    A(i, v.second->id) = v.first;
                }

                b(i,0) = constraints[i].rhs;
            }

            Matrix c(1, variables.size());
            for (const auto& v: objective_.terms) {
                if (v.first == 0) { continue; }
                c(0, v.second->id) = (type == ProblemType::Max) ? -v.first :  v.first;
            }

            // Solve
            RevisedSimplex simplex(A, b, c);
            RevisedSimplex::Solution s = simplex.solve();

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

            return os;
        }
    private:
        ProblemType type;
        Expression objective_;
        std::vector<Variable> variables;
        std::vector<Constraint> constraints;
    };
}

