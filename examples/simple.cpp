#include <iostream>
#include <vector>
#include "../lp.h"

template<typename T>
void print_vector(const std::vector<T>& v, const std::string& name = "") {
    if (!name.empty()) std::cout << name << " = ";
    for (const auto& x : v) {
        std::cout << x << " ";
    }
    std::cout << "\n";
}

using namespace lp;

size_t g_iter_n = 0;

void print_state(RevisedSimplex::State& state) {
    std::vector<lp::Number> bfs(state.A.cols(), 0);
    std::vector<lp::Number> c(state.A.cols(), 0);
    for (size_t i = 0; i < state.x.rows(); ++i) {
        bfs[i] = state.x(i, 0); 
    }

    std::cout << std::endl << "Iteration #" << g_iter_n << std::endl;
    std::cout << "-------------------" << std::endl;
    std::cout << "Status: " << state.status << std::endl;
    print_vector(state.bv, "bv");
    print_vector(state.nbv, "nbv");
    print_vector(bfs, "bfs");
    std::cout << "z = " << (state.c * state.x)(0,0) << std::endl;

    std::cout << "c: " << state.c;
    std::cout << "Binv: " << std::endl << state.Binv << std::endl;

}

void hook(RevisedSimplex::State& state) {
    g_iter_n++;
    print_state(state);
}

int main() {
    Matrix c(1, 2);
    // Obj: -x1 - x2
    c(0, 0) = -1;
    c(0, 1) = -1;

    Matrix A(2, 2);
    Matrix b(2, 1);
    // C1: 2x1 + x2 <= 4
    A(0, 0) = 2;
    A(0, 1) = 1;

    b(0, 0) = 4;

    // C2: 3x1 + 5x2 <= 15
    A(1, 0) = 3;
    A(1, 1) = 5;
    b(1, 0) = 15;

    std::cout << "c: " << c;
    std::cout << "[A | b]:" << std::endl;
    std::cout << Matrix::augment(A, b) << std::endl;

    std::cout << "Performing simplex:" << std::endl;
    g_iter_n = 0;
    RevisedSimplex splx(A, b, c);
    splx.set_hook(hook);

    RevisedSimplex::Solution soln = splx.solve();
    std::cout << std::endl << "Solution: " << std::endl;
    std::cout << "-------------------" << std::endl;
    print_vector(soln.x, "x");
    std::cout << "z* = " << soln.z << std::endl;

    return 0;
}
