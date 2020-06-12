#include <iostream>
#include "eigenvalues.h"

using namespace std;

int main()
{
    double matrix[3][3] = {
    { -1.48213, -0.03916,  1.08254 },
    { -0.03916,  1.13958,  0.01617 },
    {  1.08254,  0.01617, -1.48271 }
    };
    cout << "POWER METHOD RESULTS" << endl;
    power_method(matrix);

    cout << endl << "SCALAR PRODUCT METHOD RESULTS" << endl;
    scalar_product_method(matrix);

    return 0;
}
