#include <iostream>
#include "tridiag_matrix_method.h"
#include <cmath>
using namespace std;

double p (double x) {
    return 1.0 / (x - 3);
}

double q (double x) {
    return 1 + x / 2.0;
}

double r (double x) {
    return exp(x / 2.0);
}

double f (double x) {
    return x - 2;
}

int main()
{
    double* y10 = tridiagonal_matrix(p, q, r, f, 10, -1, 1);
    cout << "RESULT:" << endl;
    double h = 0.2;
    double x = -1;
    cout << "x\t    y10\t" << endl;
    for (int i = 0; i < 11; i++) {
        cout << x << "\t" << y10[i] << endl;
        x += h;
    }
    return 0;
}
