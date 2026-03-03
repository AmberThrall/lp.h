#include <iomanip>
#include <iostream>

#define LP_H_EIGEN
#include "../lp.h"
using namespace lp;

struct Vitamin {
    std::string name;
    double min;
    double max;
};

const Vitamin VITAMINS[] = {
    { "A", 700, 10000 },
    { "C", 700, 10000 },
    { "B1", 700, 10000 },
    { "B2", 700, 10000 },
};

struct Food {
    std::string name;
    double cost;
    double min;
    double max;
    double nutitional_value[4];
};

const Food FOODS[] = {
    { "beef",       3.19, 0, 100, { 60, 20, 10, 15} },
    { "chicken",    2.59, 0, 100, {  8,  0, 20, 20} },
    { "fish",       2.29, 0, 100, {  8, 10, 15, 10} },
    { "ham",        2.89, 0, 100, { 40, 40, 35, 10} },
    { "mac&cheese", 1.89, 0, 100, { 15, 35, 15, 15} },
    { "meatloaf",   1.99, 0, 100, { 70, 30, 15, 15} },
    { "spaghetti",  1.99, 0, 100, { 25, 50, 25, 15} },
    { "turkey",     2.49, 0, 100, { 60, 20, 15, 10} },
};

int main() {
    // Construct the problem
    Problem problem = Problem::minimize();

    std::vector<Variable> vars;
    for (auto & food : FOODS) {
        vars.push_back(problem.add_var(food.min, food.max, food.name));
    }

    Expression objective;
    for (size_t i = 0; i < vars.size(); ++i) {
        objective += FOODS[i].cost * vars[i];
    }
    problem.objective(objective);

    for (size_t i = 0; i < 4; ++i) {
        Expression expr;
        for (size_t j = 0; j < vars.size(); ++j) {
            expr += FOODS[j].nutitional_value[i] * vars[j];
        }
        problem.constraint(expr >= VITAMINS[i].min);
        problem.constraint(expr <= VITAMINS[i].max);
    }
    
    std::cout << "Problem: " << std::endl << problem;

    Solution soln = problem.solve();
    std::cout << std::endl << "Solution: ";
    if (soln.status != SolutionStatus::kOptimal) {
        std::cout << soln.status << std::endl;
    }
    else {
        std::cout << std::fixed << std::setprecision(2);
        std::cout << std::endl;
        std::cout << "Food Item  Count  Price" << std::endl;
        std::cout << "---------  -----  -----" << std::endl;
        for (size_t i = 0; i < vars.size(); ++i) {
            std::cout << FOODS[i].name;

            size_t spaces = 11 - FOODS[i].name.length();
            for (size_t j = 0; j < spaces; ++j) { std::cout << " "; }

            std::cout << std::setw(5) << soln[vars[i]];
            std::cout << "  " << std::setw(5) << soln[vars[i]] * FOODS[i].cost << std::endl;
        }
        std::cout << std::endl << "Total cost: " << soln.objective << std::endl;
    }
}
