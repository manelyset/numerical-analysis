#include <iostream>
#include <cmath>
#include "precise.h"
using namespace std;

const double accuracy = 0.0001;

/*double abs (double x) {
    return x >= 0 ? x : -x;
}*/

void multiply_vector_by_matrix(double matrix[3][3], double vect[3]) {
    double* result = new double[3];
    for (int i = 0; i < 3; i++)
        result[i] = 0;
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++) {
            result[i] += matrix[i][j] * vect[j];
        }
    for (int i = 0; i < 3; i++)
        vect[i] = result[i];
    delete[] result;

}
void add_vector(double v1[3], double v2[3]) {
    for (int i = 0; i < 3; i++)
        v1[i] += v2[i];
}
void sub_vector(double v1[3], double v2[3]) {
    for (int i = 0; i < 3; i++)
        v1[i] -= v2[i];
}
void copy_vector(double v1[3], double v2[3]) {
    for (int i = 0; i < 3; i++)
        v1[i] = v2[i];
}
double vector_norm (double vect[3]) {
    double sum = 0;
    for (int i = 0; i < 3; i++)
        sum += abs(vect[i]);
    return sum;
}
double matrix_norm (double matrix[3][3]) {
    double max_sum = 0, sum;
    for (int i = 0; i < 3; i++) {
        sum = 0;
        for (int j = 0; j < 3; j++) {
            sum += abs(matrix[i][j]);
        }
        if (sum > max_sum) max_sum = sum;
    }
    return max_sum;
}

void print_vector (double vect[3]) {
    for (int i = 0; i < 3; i++) {
        cout << vect[i] << " ";
    }
    cout << endl;
}
void simple_iteration_3x3(double a[3][4]) {
    double h[3][3];
    double g[3];
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            if (i == j)
                h[i][j] = 0;
            else
                h[i][j] = - a[i][j] / a[i][i];
        }
        g[i] = a[i][3] / a[i][i];
    }

    double nh = matrix_norm(h), ng = vector_norm(g);
    cout << "Norm of H: " << nh << endl;

    double iterations_est = log(nh) / log (accuracy * (1 - nh) / ng);
    cout << "A priori iterations estimation: " << floor(iterations_est) << endl;

    double x[3] = {0, 0, 0};
    double x1[3];
    double* precise = hauss_3x3(a);
    int iterations_count = 0;

    copy_vector(x1, x);
    sub_vector(x1, precise);
    while (vector_norm(x1) >= accuracy && iterations_count < 10000) {
        multiply_vector_by_matrix(h, x);
        add_vector(x, g);

        copy_vector(x1, x);
        sub_vector(x1, precise);
        iterations_count++;
    }
    cout << "\nResult of simple iteration method:\n";
    for (int i = 0; i < 3; i++) {
        cout << "x" << i + 1 << " = " << x[i] << endl;
    }
    cout << "Actual number of iterations: " << iterations_count << endl;

}

void seidel_3x3(double a[3][4]) {
    double h[3][3];
    double g[3];
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            if (i == j)
                h[i][j] = 0;
            else
                h[i][j] = - a[i][j] / a[i][i];
        }
        g[i] = a[i][3] / a[i][i];
    }
    cout << "\nNorm of H: " << matrix_norm(h) << endl;
    double x[3] = {0, 0, 0};
    double x1[3];
    double* precise = hauss_3x3(a);
    int iterations_count = 0;
    double xi;

    copy_vector(x1, x);
    sub_vector(x1, precise);
    while (vector_norm(x1) >= accuracy && iterations_count < 10000) {
        for (int i = 0; i < 3; i++) {
            xi = 0;
            for (int j = 0; j < 3; j++) {
                xi += h[i][j] * x[j];
            }
            xi += g[i];
            x[i] = xi;
        }
        copy_vector(x1, x);
        sub_vector(x1, precise);
        iterations_count++;
    }
    cout << "\nResult of Seidel method:\n";
    for (int i = 0; i < 3; i++) {
        cout << "x" << i + 1 << " = " << x[i] << endl;
    }
    cout << "Actual number of iterations: " << iterations_count << endl;

}

void upper_relaxation_3x3(double a[3][4]) {
    double h[3][3];
    double g[3];
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            if (i == j)
                h[i][j] = 0;
            else
                h[i][j] = - a[i][j] / a[i][i];
        }
        g[i] = a[i][3] / a[i][i];
    }
    double q = 1;


    double x[3] = {0, 0, 0};
    double x1[3];
    double* precise = hauss_3x3(a);
    int iterations_count = 0;
    double xi;

    copy_vector(x1, x);
    sub_vector(x1, precise);
    while (vector_norm(x1) >= accuracy && iterations_count < 10000) {
        for (int i = 0; i < 3; i++) {
            xi = 0;
            for (int j = 0; j < 3; j++) {
                xi += q * h[i][j] * x[j];
            }
            xi += q * g[i] + (1 - q) * x[i];
            x[i] = xi;
        }
        copy_vector(x1, x);
        sub_vector(x1, precise);
        iterations_count++;
    }
    cout << "\nResult of upper relaxation method:\n";
    for (int i = 0; i < 3; i++) {
        cout << "x" << i + 1 << " = " << x[i] << endl;
    }
    cout << "With q = " << q <<endl;
    cout << "Actual number of iterations: " << iterations_count << endl;

}

