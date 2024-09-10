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

#include "linearoscp_gvf.h"


int main (int argc, char *argv[])
{
    const int N = 2;
    double z[N];
    double params[1] = {1.0};

    const gsl_odeiv2_step_type *T  = gsl_odeiv2_step_rk8pd;
    gsl_odeiv2_step    *step    = gsl_odeiv2_step_alloc(T, N);
    gsl_odeiv2_control *control = gsl_odeiv2_control_y_new(1e-12, 1e-9);
    gsl_odeiv2_evolve  *evolve  = gsl_odeiv2_evolve_alloc(N);
    gsl_odeiv2_system sys = {linearoscp_vf, linearoscp_jac, N, params};

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
        /* Something went wrong while solving. */
        return -1;
    }

    /* Verify that we got the expected final values. */
    if (!((fabs(z[0] + 1) < 1e-9) && (fabs(z[1]) < 1e-9))) {
        fprintf(stderr, "FAIL: tolerance for final point not met\n");
        return 1;
    }
    double sqmag = linearoscp_func(0.0, z, params);
    if (fabs(sqmag - 1) > 1e-10) {
        fprintf(stderr, "FAIL: squared magnitude should be approx. 1.\n");
        fprintf(stderr, "      got sqmag = %25.18e\n", sqmag);
        return 1;
    }
    double theta = linearoscp_theta(0.0, z, params);
    if (fabs(theta - M_PI) > 1e-10) {
        fprintf(stderr, "FAIL: theta should be approx. pi\n");
        fprintf(stderr, "      got theta = %25.18e\n", theta);
        return 1;
    }
    double f[2];
    linearoscp_functions(0.0, z, params, f);
    if (f[linearoscp_function_func] != sqmag ||
            f[linearoscp_function_theta] != theta) {
        fprintf(stderr, "FAIL: values returned by linearoscp_functions do not match.\n");
        fprintf(stderr, "      the values returned by the individual functions.\n");
        fprintf(stderr, "      sqmag = %25.18e\n", sqmag);
        fprintf(stderr, "      f[0]  = %25.18e\n", f[0]);
        fprintf(stderr, "      theta = %25.18e\n", theta);
        fprintf(stderr, "      f[1]  = %25.18e\n", f[1]);
    }
}
