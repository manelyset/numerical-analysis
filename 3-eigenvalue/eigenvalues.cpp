#include <iostream>
#define accuracy 0.001
#include <cmath>

using namespace std;

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
double vector_norm (double vect[3]) {
    double sum = 0;
    for (int i = 0; i < 3; i++)
        sum += vect[i] * vect[i];
    return sqrt(sum);
}

void copy_vector(double v1[3], double v2[3]) {
    for (int i = 0; i < 3; i++)
        v1[i] = v2[i];
}

void normalize (double v[3]) {
    double norm = vector_norm(v);
    for (int i = 0; i < 3; i++) {
        v[i] /= norm;
    }
}

double scalar_product (double v1[3], double v2[3]) {
    double result = 0;
    for (int i = 0; i < 3; i++) {
        result += v1[i] * v2[i];
    }
    return result;
}

void scalar_product_method (double matrix[3][3]) {
    double y_old[3] = {1, 1, 1};
    double y_new[3] = {1, 1, 1};
    double value_old, value_new = 1;
    int iterations_number = 0;
    //normalize(y_old);
    normalize(y_new);
    do{
        value_old = value_new;
        copy_vector(y_old, y_new);
        multiply_vector_by_matrix(matrix, y_new);
        value_new = scalar_product(y_old, y_new) / scalar_product(y_old, y_old);
        normalize(y_new);
        iterations_number++;
    }
    while (fabs(value_new - value_old) > accuracy);
    cout << "Eigen value: " << value_new << endl;
    cout << "Eigen vector: [ ";
    for (int i = 0; i < 3; i++) {
        cout << y_new[i] << " ";
    }
    cout << "]" << endl;
    cout << "Number of iterations: " << iterations_number << endl;
}

void power_method (double matrix[3][3]) {
    double y_old[3] = {1, 1, 1};
    double y_new[3] = {1, 1, 1};
    double value;
    int iterations_number = 0;
    normalize(y_new);
    do {
        copy_vector(y_old, y_new);
        multiply_vector_by_matrix(matrix, y_new);
        value = y_new[0] / y_old[0];
        normalize(y_new);
        iterations_number++;
    } while (fabs(y_old[0] / y_new[0] - y_old[1] / y_new[1]) > accuracy || fabs(y_old[2] / y_new[2] - y_old[1] / y_new[1]) > accuracy);
    cout << "Eigen value: " << value << endl;
    cout << "Eigen vector: [ ";
    for (int i = 0; i < 3; i++) {
        cout << y_new[i] << " ";
    }
    cout << "]" << endl;
    cout << "Number of iterations: " << iterations_number << endl;
}
