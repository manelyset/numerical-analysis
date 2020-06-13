#include "fredholm.h"
#include <cmath>
#include <iostream>
#include <vector>
#include "slau_solving.h"
#define accuracy 0.001

using namespace std;

double f (double x) {
    return x - 0.6;
}

double H (double x, double y) {
    return 0.6 * cosh(x * y);
}

int factorial (int x) {
    return x <= 1 ? 1 : x * factorial(x - 1);
}



double sum (double x, double* c, int n) {
    double result = 0;
    for (int i = 0; i < n; i++) {
        result += c[i] * (pow(x, 2.0 * i) / factorial(2 * i));
    }
    return result;
}

void degenerate_kernel() {
    double gamma[4][4], b[4];

    for (int i = 0; i < 4; i++) {
        b[i] = 1 / (2 * i + 2) - 0.6 / (2 * i + 1);
        for (int j = 0; j < 4; j++) {
            gamma[i][j] = 1.0 / (factorial(2*j) * (i * j + 1));
        }
    }

    double a3[3][4];
    double a4[4][5];

    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            a3[i][j] = -gamma[i][j];
        }
        a3[i][i] += 1;
        a3[i][3] = b[i];
    }
    swap_lines(a3, 0, 1);

    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            a4[i][j] = -gamma[i][j];
        }
        a4[i][i] += 1;
        a4[i][4] = b[i];
    }
    swap_lines_4(a4, 0, 1);

    double* c3 = jordan_3x3(a3);
    double* c4 = jordan_4x4(a4);

    double u3[3] = { f(0) + sum (0, c3, 3) , f(1/2) + sum(1/2, c3, 3) , f(1) + sum(1, c3, 3) };
    double u4[3] = { f(0) + sum (0, c4, 4) , f(1/2) + sum(1/2, c4, 4) , f(1) + sum(1, c4, 4) };



    cout << "DEGENERATE KERMEL METHOD" << endl;
    cout << "func\t0\t1/2\t1" << endl;
    cout << "u_3\t";
    for (int i = 0; i < 3; i++)
        cout << u3[i] << "\t";
    cout << endl << "u_4\t";
     for (int i = 0; i < 3; i++)
        cout << u4[i] << "\t";
    cout << endl;
    double max_diff = 0;
    for (int i = 0; i < 3; i++) {
        if (fabs(u3[i] - u4[i]) > max_diff)
            max_diff = fabs(u3[i] - u4[i]);
    }
    cout << "Estimated max difference: " << max_diff << endl;

}





void mechanic_quadrature() {
    int n = 2;
    double Ak = 1.0 / n;

    double result_old[3], result_new[3] = {0, 0, 0};

    do {
        vector<vector<double> > D(n, vector<double>(n + 1));
        copy_array(result_old, result_new);
        //cout << "u\t" << result_new[0] << "\t" << result_new[1] <<"\t" << result_new[2] <<endl;

        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                D[i][j] = -Ak * H((double)i/n + (double) 1/(2 * n), (double)j/n + (double) 1/(2 * n));
                }
            D[i][i] += 1;
            D[i][n] = f((double)i/n + (double) 1/(2 * n));

        }

        vector <double> z = jordan(D, n);
        result_new[0] = f(0);
        for (int i = 0; i < n; i++) {
            result_new[0] += Ak * H(0, (double)i/n + (double) 1/(2 * n)) * z[i];
        }
        result_new[1] = f(1/2.0);
        for (int i = 0; i < n; i++) {
            result_new[1] += Ak * H(1/2.0, (double)i/n + (double) 1/(2 * n)) * z[i];
        }
        result_new[2] = f(1);
        for (int i = 0; i < n; i++) {
            result_new[2] += Ak * H(1, (double)i/n + (double) 1/(2 * n)) * z[i];
        }
        D.clear();
        z.clear();
        n++;
        Ak = 1.0/n;
        //cout << n << endl;
    } while ((fabs(result_new[0] - result_old[0]) > accuracy || fabs(result_new[1] - result_old[1]) > accuracy || fabs(result_new[2] - result_old[2]) > accuracy)&& n < 50);


    cout << "MECHANIC QUADRATURE METHOD" << endl;
    cout << "x\t0\t1/2\t1" << endl;
    cout << "u\t" << result_new[0] << "\t" << result_new[1] <<"\t" << result_new[2] <<endl;

}
