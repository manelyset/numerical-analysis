#ifndef FINITE_DIFFERENCES_H_INCLUDED
#define FINITE_DIFFERENCES_H_INCLUDED

void explicit_method (int n, int m, double (*f)(double, double), double T, double (*phi)(double), double (*alpha)(double), double (*beta)(double));
void implicit_method (int n, int m, double (*f)(double, double), double T, double (*phi)(double), double (*alpha)(double), double (*beta)(double));

#endif // FINITE_DIFFERENCES_H_INCLUDED
