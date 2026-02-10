#include <iostream>

#include "../lp.h"
using namespace lp;

int main() {
    Problem problem = Problem::maximize();
    Variable x1 = problem.add_var(0, 5);
    Variable x2 = problem.add_var(-10, 10);
    problem.objective(30*x1 - 4*x2);
    problem.constraint(5*x1 - x2 <= 30);

    std::cout << "Problem: " << std::endl << problem;

    Solution soln = problem.solve();
    std::cout << std::endl << "Solution: ";
    if (soln.status != SolutionStatus::kOptimal) {
        std::cout << soln.status << std::endl;
    }
    else {
        std::cout << std::endl;
        std::cout << "x1  = " << soln[x1] << std::endl;
        std::cout << "x2  = " << soln[x2] << std::endl;
        std::cout << "z* = " << soln.objective << std::endl;
    }

    return 0;
}
