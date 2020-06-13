#include "fredholm.h"
#include <cmath>
#include <iostream>
#include<vector>
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

void swap_lines_4 (double matrix[4][5], int first, int second) {
    for (int i = 0; i <= 4; i++)
        swap (matrix[first][i], matrix[second][i]);
}
void divide_line_4 (double matrix[4][5], int j, double x) {
    for (int i = 0; i <= 4; i++)
        matrix[j][i] /= x;
}
void sub_lines_4 (double matrix[4][5], int a, int b, double k) {
    for (int i = 0; i <= 4; i++)
        matrix[a][i] -= matrix[b][i] * k;
}

double* jordan_4x4 (double matrix1[4][5]) {
    double matrix[4][5];
    for (int i = 0; i < 4; i++)
        for (int j = 0; j < 5; j++)
            matrix[i][j] = matrix1[i][j];
    for (int i = 0; i < 4; i++) {
        divide_line_4(matrix, i, matrix[i][i]);
    }

    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            if (i != j) {
                sub_lines_4(matrix, j, i, matrix[j][i]);
            }
        }
    }
    double* answer = new double[4];
    for (int i = 0; i < 4; i++) {
        answer[i] = matrix[i][4];
    }
    return answer;
}
void swap_lines (double matrix[3][4], int first, int second) {
    for (int i = 0; i <= 3; i++)
        swap (matrix[first][i], matrix[second][i]);
}
void divide_line (double matrix[3][4], int j, double x) {
    for (int i = 0; i <= 3; i++)
        matrix[j][i] /= x;
}
void sub_lines (double matrix[3][4], int a, int b, double k) {
    for (int i = 0; i <= 3; i++)
        matrix[a][i] -= matrix[b][i] * k;
}
double* jordan_3x3 (double matrix1[3][4]) {
    double matrix[3][4];
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 4; j++)
            matrix[i][j] = matrix1[i][j];
    for (int i = 0; i < 3; i++) {
        divide_line(matrix, i, matrix[i][i]);
    }

    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            if (i != j) {
                sub_lines(matrix, j, i, matrix[j][i]);
            }
        }
    }
    double* answer = new double[3];
    for (int i = 0; i < 3; i++) {
        answer[i] = matrix[i][3];
    }
    return answer;
}

double sum (double x, double* c, int n) {
    double result = 0;
    for (int i = 0; i < n; i++) {
        result += c[i] * (pow(x, 2.0 * i) / factorial(2 * i));
    }
    return result;
}


void degenerate_kernel() {
    /*ch(xy) = 1*1 + (xy)^2/2! + (xy)^4/4! + (xy)^6/6!
    alpha_i(x) = x^(2i)/(2i)!
    beta_i(y) = y^(2i)
    */

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

void swap_lines_vector (vector<vector<double> > matrix, int first, int second, int n) {
    for (int i = 0; i <= n; i++)
        swap (matrix[first][i], matrix[second][i]);
}

void divide_line_vector (vector<vector<double> > matrix, int j, double x, int n) {
    for (int i = 0; i <= n; i++)
        matrix[j][i] /= x;
}
void sub_lines_vector (vector<vector<double> > matrix, int a, int b, double k , int n) {
    for (int i = 0; i <= n; i++)
        matrix[a][i] -= matrix[b][i] * k;
}
void copy_array(double a1[3], double a2[3]) {
    for (int i = 0; i < 3; i++)
        a2[i] = a1[i]; }

vector<double> hauss (vector<vector<double> > matrix, int n) {

    for (int i = 0; i < n; i++) {
        if (matrix[i][i] == 0) {
            for (int j = i; j < n; j++) {
                if (matrix[j][i] != 0) swap_lines_vector(matrix, i, j, n);
            }
        }
        divide_line_vector(matrix, i, matrix[i][i], n);
        for (int j = i; j < n; j++)
            sub_lines_vector(matrix, j, i, matrix[i][j], n);
    }


    vector<double> answer(n);
    for (int i = n - 1; i >= 0; i--) {
        answer[i] = matrix[i][n+1];
        for (int j = i + 1; j < n; j++)
            answer[i] -= matrix[i][j] * answer[j];
    }
    return answer;

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

        vector <double> z = hauss(D, n);
        result_new[0] = f(0);
        for (int i = 0; i < n; i++) {
            result_new[0] += Ak * H(0, (double)i/n + (double) 1/(2 * n)) * z[i];
        }
        result_new[1] = f(1/2);
        for (int i = 0; i < n; i++) {
            result_new[1] += Ak * H(1/2.0, (double)i/n + (double) 1/(2 * n)) * z[i];
        }
        result_new[2] = f(1);
        for (int i = 0; i < n; i++) {
            result_new[2] += Ak * H(1, (double)i/n + (double) 1/(2 * n)) * z[i];
        }
        n *= 2;
        Ak = 1.0/n;
        //cout << n << endl;
    } while ((fabs(result_new[0] - result_old[0]) > accuracy || fabs(result_new[1] - result_old[1]) > accuracy || fabs(result_new[2] - result_old[2]) > accuracy)&& n < 128);


    cout << "MECHANIC QUADRATURE METHOD" << endl;
    cout << "x\t0\t1/2\t1" << endl;
    cout << "u\t" << result_new[0] << "\t" << result_new[1] <<"\t" << result_new[2] <<endl;

}
