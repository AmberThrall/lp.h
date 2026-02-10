#include <iostream>

#include "../lp.h"
using namespace lp;

class TraceSolver : public DefaultSolver {
public:
    Solution solve(Matrix A, Matrix b, Matrix c) override {
        std::cout << "Performing revised simplex method:" << std::endl;

        return DefaultSolver::solve(A, b, c); 
    };
protected:
    void start() override {
        DefaultSolver::start();
        print_status();
    }

    void step() override {
        DefaultSolver::step();
        print_status();
    }
private:
    void print_status() {
        std::cout << "Iteration #" << iter_num << ": ";
        std::cout << "x = <";
        for (size_t i = 0; i < x.rows(); ++i) { 
            if (i > 0) { std::cout << ","; }
            std::cout << x(i,0);
        }
        std::cout << ">; z* = " << obj_value() << " (" << status << ")" << std::endl;

    }
};

int main() {
    Problem problem = Problem::minimize();
    Variable x = problem.add_var(0, lp::Infinity, "x");
    Variable y = problem.add_var(0, lp::Infinity, "y");
    problem.objective(-x - y);
    problem.constraint(2*x + y <= 4);
    problem.constraint(3*x + 5*y <= 15);

    std::cout << "Problem: " << std::endl << problem;
    std::cout << std::endl;

    TraceSolver * solver = new TraceSolver();
    Solution soln = problem.solve(solver);
    delete solver;

    std::cout << std::endl << "Solution: ";
    if (soln.status != SolutionStatus::kOptimal) {
        std::cout << soln.status << std::endl;
    }
    else {
        std::cout << std::endl;
        std::cout << "x  = " << soln[x] << std::endl;
        std::cout << "y  = " << soln[y] << std::endl;
        std::cout << "z* = " << soln.objective << std::endl;
    }

    std::cout << std::endl;
}
