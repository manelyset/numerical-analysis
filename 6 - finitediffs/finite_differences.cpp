#include "finite_differences.h"
#include <iostream>
using namespace std;

void explicit_method (int n, int m, double (*f)(double, double), double T, double (*phi)(double), double (*alpha)(double), double (*beta)(double)) {
    double h = 1.0 / n;
    double tau = T / m;
    double** u = new double* [m + 1];
    for (int i = 0; i <= m; i++) {
        u[i] = new double [n + 1];
    }
    for (int i = 0; i <= n; i++) {
        u[i][0] = phi(h * i);
    }

    for (int k = 1; k <= m; k++){
        for (int i = 1; i < n; i++) {
            u[i][k] = (u[i][k - 1] + u[i + 1][k - 1] - 2 * u[i][k - 1] + u[i - 1][k - 1]) / (h * h) + (u[i + 1][k - 1] - u[i - 1][k - 1]) / (2 * h) + f(h * i, tau * (k - 1));
        }
        u[0][k] = alpha(tau * k);
        u[n][k] = beta(tau * k) * 2 * h + 4 * u[n - 1][k] - u[n - 2][k];
    }
    cout << "x\\t\t";
    for (int i = 0; i <= n; i++) {
        cout << i * h << "\t";
    }
    cout << endl << "---------------------------------------------------------------------------------"<< endl;
    for (int k = 0; k <= m; k++) {
        cout << k * tau << "\t|\t";
        for (int i = 0; i <= n; i++)
            cout << u[i][k] << "\t";
        cout << endl;
    }
}

void implicit_method (int n, int m, double (*f)(double, double), double T, double (*phi)(double), double (*alpha)(double), double (*beta)(double)) {
    double h = 1.0 / n;
    double tau = T / m;
    double sigma = 0.5;
    double** u = new double* [m + 1];
    for (int i = 0; i <= m; i++) {
        u[i] = new double [n + 1];
    }
    for (int i = 0; i <= n; i++) {
        u[i][0] = phi(h * i);
    }
    double* s = new double[n + 1];
    double* t = new double[n + 1];
    for (int k = 2; k <= m; k++) {
        double C = sigma/(h*h) + sigma/(2*h);
        double B = - 2 * sigma / (h * h) - 1 / tau;
        double A = sigma/(h*h) - sigma/(2*h);
        double G = -u[0][k-1] / tau -  + f(0, tau * k);
        s[0] = C / B;
        t[0] = -G / B;
        for (int i = 1; i < n; i++) {
            G = -u[i][k-1] / tau -  + f(i * h, tau * k);
            s[i] = C / (B - A * s[i-1]);
            t[i] = (A * t[i - 1] - G) / (B - A * s[i-1]);
        }
        u[n][k] = t[n];
        for (int i = n - 1; i >= 1; i--) {
            u[i][k] = s[i] * u[i + 1][k] + t[i];
        }
    }
    cout << "x\\t\t";
    for (int i = 0; i <= n; i++) {
        cout << i * h << "\t";
    }
    cout << endl << "---------------------------------------------------------------------------------"<< endl;
    for (int k = 0; k <= m; k++) {
        cout << k * tau << "\t|\t";
        for (int i = 0; i <= n; i++)
            cout << u[i][k] << "\t";
        cout << endl;
    }


}


