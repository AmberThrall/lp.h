# lp.h

lp.h is a drop-in C++, single-file header-only library for linear programming.

## Features
- Easy to use interface for constructing LP problems.
- Automatic conversion to standard form.
- Built-in revised simplex method solver.
- Framework for using alternative solvers.

## Installation

Simply copy `lp.h` into your project's directory and include it.

## Example

```c++
#include <iostream>
#include "lp.h"
using namespace lp;

int main() {
    Problem problem = Problem::minimize();
    Variable x = problem.add_var(0, lp::Infinity, "x");
    Variable y = problem.add_var(0, lp::Infinity, "y");
    problem.objective(-x - y);
    problem.constraint(2*x + y <= 4);
    problem.constraint(3*x + 5*y <= 15);

    std::cout << "Problem: " << std::endl << problem;
    std::cout << std::endl;

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
}
```
