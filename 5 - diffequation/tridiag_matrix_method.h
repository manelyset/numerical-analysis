#ifndef TRIDIAG_MATRIX_METHOD_H_INCLUDED
#define TRIDIAG_MATRIX_METHOD_H_INCLUDED

double* tridiagonal_matrix (double (*p)(double), double (*q)(double), double (*r)(double), double (*f)(double), int n, double a, double b);

#endif // TRIDIAG_MATRIX_METHOD_H_INCLUDED
