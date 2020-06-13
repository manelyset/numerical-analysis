#include "slau_solving.h"
#include <vector>

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

vector<double> jordan (vector<vector<double> > matrix, int n) {

    for (int i = 0; i < n; i++) {
        divide_line_vector(matrix, i, matrix[i][i], n);
    }

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            if (i != j) {
                sub_lines_vector(matrix, j, i, matrix[j][i], n);
            }
        }
    }
    vector <double> answer(n);
    for (int i = 0; i < n; i++) {
        answer[i] = matrix[i][n];
    }
    return answer;
}
