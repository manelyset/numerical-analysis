#include <iostream>
#include "finite_differences.h"

using namespace std;

double f1 (double x, double t) {
    return 3 * t * t + 6 * x - 3 * x * x;
}

double phi1 (double x) {
    return x * x * x;
}

double alpha1 (double t) {
    return t * t * t;
}

double beta1 (double t) {
    return 0;
}

int main()
{
    cout << "EXPLICIT METHOD" << endl;
    explicit_method(5, 5, f1, 0.1, phi1, alpha1, beta1);
    cout << endl << "IMPLICIT METHOD, sigma = 1" << endl;
    implicit_method(5, 5, f1, 0.1, phi1, alpha1, beta1);
    return 0;
}
