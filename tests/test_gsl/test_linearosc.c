/*
 *  test_gsl.c
 *
 *
 *  GSL ODE solver for the vector field named: linearosc
 *
 *  To compile and run this program:
 *      gcc -c linearosc_gvf.c
 *      gcc -c test_gsl.c
 *      gcc -o test_gsl test_gsl.o linearosc_gvf.o -lgsl -lgslcblas -lm
 *  This creates an executable file called test_gsl.
 */

#include <string.h>
#include <math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv2.h>

#include "linearosc_gvf.h"


int main (int argc, char *argv[])
{
    void *params = NULL;
    const int N = 2;
    double z[N];

    const gsl_odeiv2_step_type *T  = gsl_odeiv2_step_rk8pd;
    gsl_odeiv2_step    *step    = gsl_odeiv2_step_alloc(T, N);
    gsl_odeiv2_control *control = gsl_odeiv2_control_y_new(1e-12, 1e-9);
    gsl_odeiv2_evolve  *evolve  = gsl_odeiv2_evolve_alloc(N);
    gsl_odeiv2_system sys = {linearosc_vf, linearosc_jac, N, params};

    double t  = 0.0;
    double t1 = M_PI;
    double h = 1e-6;
    z[0] = 1.0;
    z[1] = 0.0;
    int fail = 0;
    while (t < t1) {
        int status = gsl_odeiv2_evolve_apply(evolve, control, step,
                                             &sys, &t, t1, &h, z);
        if (status != GSL_SUCCESS) {
            fprintf(stderr, "status=%d\n", status);
            fail = 1;
            break;
        }
    }

    gsl_odeiv2_evolve_free(evolve);
    gsl_odeiv2_control_free(control);
    gsl_odeiv2_step_free(step);

    if (fail) {
        return -1;
    }
    if ((fabs(z[0] + 1) < 1e-9) && (fabs(z[1]) < 1e-9)) {
        return 0;
    }
    else {
        fprintf(stderr, "FAIL: tolerance for final point not met\n");
        return 1;
    }
}
