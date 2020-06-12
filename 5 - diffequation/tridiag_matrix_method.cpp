#include "tridiag_matrix_method.h"
#include <cmath>
#include <iostream>

double* tridiagonal_matrix (double (*p)(double), double (*q)(double), double (*r)(double), double (*f)(double), int n, double a, double b) {
    double h = (b - a) / n; //отрезок (-1, 1)
    double* s = new double[n + 1];
    double* t = new double[n + 1];
    s[0] = 0;
    t[0] = 0;
    for (int i = 1; i < n; i++) {
        double node = i * h + a;
        s[i] = ((-p(node))/(h * h) - q(node) / (2 * h)) / (((2 * p(node))/(h * h) + r(node)) - ((-p(node))/(h * h) + q(node) / (2 * h)) * s[i - 1]);
        t[i] = (((-p(node))/(h * h) + q(node) / (2 * h)) * t[i - 1] - f(node)) / (((2 * p(node))/(h * h) + r(node)) - ((-p(node))/(h * h) + q(node) / (2 * h)) * s[i - 1]);
    }
    s[n] = 0;
    t[n] = 0;
    double* y = new double[n + 1];
    y[n] = 0;
    for (int i = n - 1; i > 0; i--) {
        y[i] = s[i] * y[i + 1] + t[i];
        //std::cout << i << " " << y[i] << std::endl;
    }
    return y;
}
