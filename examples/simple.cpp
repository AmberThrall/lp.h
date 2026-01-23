#include <iostream>
#include "../lp.h"

using namespace lp;

int main() {

    Matrix m(4,6);
    m(0,0) = 10;
    m(0,1) = 20;
    m(1,1) = 30;
    m(1,3) = 40;
    m(2,2) = 50;
    m(2,3) = 60;
    m(2,3) = 70;
    m(3,3) = 10;
    m(3,5) = 80;

    Matrix b = m.submatrix({0, 1, 2, 3});
    Matrix b_inv = b.inverse();


    std::cout << m << std::endl;
    std::cout << b << std::endl;
    std::cout << b_inv << std::endl;
    std::cout << b * b_inv << std::endl;

    return 0;
}
