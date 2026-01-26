#include <iostream>

#include "../lp.h"
using namespace lp;

int main() {
    Problem problem = Problem::minimize();
    auto x = problem.add_var(0, lp::Inf, "x");
    auto y = problem.add_var(0, lp::Inf, "y");
    problem.objective(-x - y);
    problem.constraint(2*x + y <= 4);
    problem.constraint(3*x + 5*y <= 15);

    std::cout << "Problem: " << std::endl << problem;

    Solution soln = problem.solve();
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

    return 0;
}
