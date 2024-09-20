/*
 *  test_linearoscp.cpp
 */

#include <iostream>
#include <iomanip>
#include <math.h>

#include "adolc/drivers/odedrivers.h"
#include "adolc/adalloc.h"

#include "linearoscp_adolc.h"

using namespace std;

#define PI 3.1415926535897932384626433

void print_solution(double t, double x[])
{
    cout << setprecision(5) << fixed << setw(8) << t << ", ";
    cout << setprecision(10) << setw(15) << x[0] << ", " << setw(15) << x[1] << endl;
}

int main (int argc, char *argv[])
{
    double p[1];
    double x[2], xnew[2];
    double f[2];
    double t, tfinal;
    double h;
    int numsteps;
    int deg = 7;
    const int tag = 1;

    double **taylorcoeffs = myalloc2(2, deg+1);  // myalloc2 is from adol-c.

    // omega
    p[0] = 0.5;

    //
    // Initial conditions.
    //
    x[0] = 1.0;   // theta(0)
    x[1] = 0.0;   // theta'(0)


    tfinal = PI;
    numsteps = 400;
    h = tfinal/numsteps;        // step size
    t = 0.0;
    linearoscp_vf(tag, x, f, p); // Call linearoscp_vf once to compute the ADOL-C data
    for (int i = 0; i < numsteps; ++i) {
        // Call forode to compute the Taylor coefficients at the current x
        taylorcoeffs[0][0] = x[0];
        taylorcoeffs[1][0] = x[1];
        forode(tag, 2, deg, taylorcoeffs);
        // Use the Taylor coefficients to compute the approximate value xnew = x(t+h)
        double hj = 1.0;
        xnew[0] = 0.0;
        xnew[1] = 0.0;
        for (int j = 0; j < deg+1; ++j) {
            for (int i = 0; i < 2; ++i) {
                xnew[i] += taylorcoeffs[i][j]*hj;
            }
            hj = hj*h;
        }
        x[0] = xnew[0];
        x[1] = xnew[1];
        t = t + h;
    }

    //
    // Final state of x should be [1, 0].
    // Verify that we got the expected final values.
    //
    if (!((fabs(x[0]) < 1e-9) && (fabs(x[1] - 1) < 1e-9))) {
        fprintf(stderr, "FAIL: tolerance for final point not met\n");
        return 1;
    }
    return 0;
}
