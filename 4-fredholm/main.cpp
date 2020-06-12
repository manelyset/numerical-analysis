#include <iostream>
#include <cmath>
#include "fredholm.h"

using namespace std;



int main()
{
    degenerate_kernel();
    cout << endl;
    mechanic_quadrature();
    return 0;
}
