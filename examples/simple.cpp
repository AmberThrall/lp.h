#include <iostream>

#define LP_H_DEBUG
#include "../lp.h"
using namespace lp;

int main() {
    Problem problem = Problem::minimize();
    Variable x = problem.add_var(0, lp::Inf, "x");
    Variable y = problem.add_var(0, lp::Inf, "y");
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

    std::cout << std::endl;

    Problem problem2 = Problem::maximize();
    x = problem2.add_var(0, lp::Inf, "x");
    y = problem2.add_var(0, lp::Inf, "y");
    problem2.objective(40*x - 30*y);
    problem2.constraint(x - y <= 12);
    problem2.constraint(2*x - y <= 16);

    std::cout << "Problem: " << std::endl << problem2;

    soln = problem2.solve();
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
