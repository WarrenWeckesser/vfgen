/*
 *  test_evf.c
 *
 */

#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdbool.h>
/*
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv.h>
*/

#include "rossler_evf_gvf.h"


int main (int argc, char *argv[])
{
    double a, b, c;
    double x, y, z;
    double dx, dy, dz;

    double params[3];
    double y0[6];
    double f[6];
    double t;
    int status, result;
    
    a = 1.0;
    b = 2.0;
    c = -1.0;

    t = 0.0;

    x = 1.0;
    y = 2.0;
    z = -2.0;
    dx = 5.0;
    dy = 2.0;
    dz = -3.0;
    
    y0[0] = x;
    y0[1] = y;
    y0[2] = z;
    y0[3] = dx;
    y0[4] = dy;
    y0[5] = dz;
    
    params[0] = a;
    params[1] = b;
    params[2] = c;

    status = rossler_evf_vf(t, y0, f, (void *)params);
    
    bool testvf = (f[0] == -(y + z)) &&
                  (f[1] == x + a*y) &&
                  (f[2] == b + (x - c)*z) &&
                  (f[3] == -dy - dz) &&
                  (f[4] == dx + a*dy) &&
                  (f[5] == z*dx + (x - c)*dz);

    if (!testvf) {
        int i;
        for (i = 0; i < 6; ++i) {
            printf("f[%d] = %20.15f\n", i, f[i]);
        }
    }

    return !testvf;
}
