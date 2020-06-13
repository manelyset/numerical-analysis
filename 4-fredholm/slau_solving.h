#ifndef SLAU_SOLVING_H_INCLUDED
#define SLAU_SOLVING_H_INCLUDED

#include <vector>
using namespace std;

void swap_lines_4 (double matrix[4][5], int first, int second);
void divide_line_4 (double matrix[4][5], int j, double x);
void sub_lines_4 (double matrix[4][5], int a, int b, double k);
double* jordan_4x4 (double matrix1[4][5]);

void swap_lines (double matrix[3][4], int first, int second);
void divide_line (double matrix[3][4], int j, double x);
void sub_lines (double matrix[3][4], int a, int b, double k);
double* jordan_3x3 (double matrix1[3][4]);

void swap_lines_vector (vector<vector<double> > matrix, int first, int second, int n);
void divide_line_vector (vector<vector<double> > matrix, int j, double x, int n);
void sub_lines_vector (vector<vector<double> > matrix, int a, int b, double k , int n);
void copy_array(double a1[3], double a2[3]);
vector<double> jordan (vector<vector<double> > matrix, int n);


#endif // SLAU_SOLVING_H_INCLUDED
