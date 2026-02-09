#include <iostream>

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
    problem2.objective(40*x + 30*y);
    problem2.constraint(x + y <= 12);
    problem2.constraint(2*x + y <= 16);

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

    // min -x - y
    // s.t. 2x + y + s1 = 4
    //      3x + 5y + s2 = 15
    /*Matrix A(2, 4);
    A(0,0) = 2;
    A(0,1) = 1;
    A(0,2) = 1;
    A(1,0) = 3;
    A(1,1) = 5;
    A(1,3) = 1;

    Matrix b(2, 1);
    b(0,0) = 4;
    b(1,0) = 15;

    Matrix c(1, 4);
    c(0,0) = -1;
    c(0,1) = -1;

    std::cout << "A =" << std::endl << A;
    std::cout << "b =" << std::endl << b;
    std::cout << "c =" << std::endl << c;

    Solver solver(A, b, c);
    Solver::Solution soln = solver.solve();
    std::cout << "soln = [";
    for (size_t i = 0; i < soln.x.size(); ++i) {
        if (i > 0) { std::cout << ","; }
        std::cout << soln.x[i];
    }
    std::cout << "], z* = " << soln.z << std::endl;*/
}
