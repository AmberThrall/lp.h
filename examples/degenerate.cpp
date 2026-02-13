#include <iostream>

#include "../lp.h"
using namespace lp;

int main() {
    Problem problem = Problem::minimize();
    Variable x = problem.add_var(0, lp::Infinity, "x");
    Variable y = problem.add_var(0, lp::Infinity, "y");
    problem.objective(2*x +y);
    problem.constraint(x + y >= 2);
    problem.constraint(3*x + y >= 4);
    problem.constraint(3*x + 2*y <= 10);
    problem.constraint(x <= 10/3.0);
    problem.constraint(y <= 5);

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
