#include <iostream>
using namespace std;
double abs (double x) {
    return x >= 0 ? x : -x;
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

double* hauss_3x3 (double matrix1[3][4]) {
    double matrix[3][4];
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 4; j++)
            matrix[i][j] = matrix1[i][j];
    if (abs(matrix[0][0]) < abs(matrix[1][0]) || abs(matrix[0][0]) < abs(matrix[2][0])) {
        if (abs(matrix[1][0]) > abs(matrix[2][0])) {
            swap_lines(matrix, 0, 1);
        }
        else {
            swap_lines(matrix, 0, 2);
        }
    }

    divide_line(matrix, 0, matrix[0][0]);
    sub_lines (matrix, 1, 0, matrix[1][0]);
    sub_lines (matrix, 2, 0, matrix[2][0]);

    if (abs(matrix[2][1]) > abs(matrix[1][1]))
        swap_lines(matrix, 1, 2);

    divide_line(matrix, 1, matrix[1][1]);
    sub_lines(matrix, 2, 1, matrix[2][1]);


    divide_line(matrix, 2, matrix[2][2]);

    double* answer = new double[3];
    answer[2] = matrix[2][3];
    answer[1] = matrix[1][3] - matrix[1][2] * answer[2];
    answer[0] = matrix[0][3] - matrix[0][2] * answer[2] - matrix[0][1] * answer[1];
    return answer;

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

double* lu_3x3 (double matrix1[3][4]) {
    double matrix[3][4];
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 4; j++)
            matrix[i][j] = matrix1[i][j];
    double l[3][4], u[3][4];
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 4; j++) {
            l[i][j] = 0;
            u[i][j] = 0;
        }
    }
    double s = 0;
    for (int i = 0; i < 3; i++) {
        for (int j = i; j < 4; j++) {
            if (j < 3) {
                for (int k = 0; k <= i - 1; k++)
                    s += l[j][k] * u[k][i];
                l[j][i] = matrix[j][i] - s;
            }
            s = 0;

            for (int k = 0; k <= i - 1; k++)
                s += l[i][k] * u[k][j];
            u[i][j] = (matrix[i][j] - s) / l[i][i];
        }
    }
    double* answer = hauss_3x3(u);
    return answer;
}
