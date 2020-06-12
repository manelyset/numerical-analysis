#include <iostream>
#include "precise.h"
#include "iterations.h"

using namespace std;

int main()
{

    double matrix[3][4] = {
    {  8.29381 , 0.995516, -0.560617, 0.766522 },
    {  0.995516, 6.298198,  0.595772, 3.844422 },
    { -0.560617, 0.595772,  4.997407, 5.239231 }
    };
    for (int i = 0; i < 3; i++) {
        cout << matrix[i][0] << " * x1 + " << matrix[i][1] << " * x2 + " << matrix[i][2] << " * x3 = " << matrix[i][3] << endl;
    }
    cout << "\n---PRECISE METHODS---\n";
    double* answer_hauss = hauss_3x3(matrix);
    cout << "\nSolved with Hauss method\n";
    for (int i = 0; i < 3; i++) {
        cout << "x" << i + 1 << " = " << answer_hauss[i] << endl;
    }
    double* answer_jordan = jordan_3x3(matrix);
    cout << "\nSolved with Jordan method\n";
    for (int i = 0; i < 3; i++) {
        cout << "x" << i + 1 << " = " << answer_jordan[i] << endl;
    }
    double* answer_lu = lu_3x3(matrix);
    cout << "\nSolved with LU method\n";
    for (int i = 0; i < 3; i++) {
        cout << "x" << i + 1 << " = " << answer_lu[i] << endl;
    }

    cout << "\n---ITERATIONS METHODS---\n";
    simple_iteration_3x3(matrix);
    seidel_3x3(matrix);
    upper_relaxation_3x3(matrix);

    return 0;
}
