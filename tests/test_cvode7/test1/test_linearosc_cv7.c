
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

/* Include headers for CVODE */
#include <cvode/cvode.h>
#include <nvector/nvector_serial.h>
#include <sunlinsol/sunlinsol_dense.h> /* dense SUNLinearSolver */
#include <sunmatrix/sunmatrix_dense.h> /* dense SUNMatrix       */

#include "linearosc_cv7.h"


int check_int_status(int retval, char *funcname)
{
    if (retval != 0) {
        fprintf(stderr, "SUNDIALS ERROR: %s() failed - returned %d\n", funcname, retval);
        return 1;
    }
    return 0;
}

int check_pointer(void *ptr, char *funcname)
{
    if (ptr == NULL) {
        fprintf(stderr, "SUNDIALS ERROR: %s() failed - returned NULL\n", funcname);
        return 1;
    }
    return 0;
}

int main (int argc, char *argv[])
{
    SUNContext sunctx;
    int i, j;
    int retval;
    const int N = 2;

    /* Create the SUNDIALS context */
    retval = SUNContext_Create(SUN_COMM_NULL, &sunctx);
    if (check_int_status(retval, "SUNContext_Create")) {return 1;}

    /* Initial conditions */
    N_Vector y = N_VNew_Serial(N, sunctx);
    if (check_pointer((void*)y, "N_VNew_Serial")) {return (1);}
    NV_Ith_S(y, 0) = SUN_RCONST(1.0);
    NV_Ith_S(y, 1) = SUN_RCONST(0.0);

    /* Use CV_ADAMS for non-stiff problems, and CV_BDF for stiff problems:   */
    void *cvode_mem = CVodeCreate(CV_ADAMS, sunctx);
    if (check_pointer((void*)cvode_mem, "CVodeCreate")) {return 1;}

    sunrealtype t = SUN_RCONST(0.0);
    retval = CVodeInit(cvode_mem, linearosc_vf, t, y);
    if (check_int_status(retval, "CVodeInit")) {return 1;}

    retval = CVodeSStolerances(cvode_mem, 5e-15, 1e-14);
    if (check_int_status(retval, "CVodeSStolerances")) {return 1;}

    /* Create dense SUNMatrix for use in linear solves */
    SUNMatrix A = SUNDenseMatrix(N, N, sunctx);
    if (check_pointer((void*)A, "SUNDenseMatrix()")) {return 1;}

    /* Create dense SUNLinearSolver object for use by CVode */
    SUNLinearSolver LS = SUNLinSol_Dense(y, A, sunctx);
    if (check_pointer((void*)LS, "SUNLinSol_Dense()")) {return 1;}

    /* Attach the matrix and linear solver */
    retval = CVodeSetLinearSolver(cvode_mem, LS, A);
    if (check_int_status(retval, "CVodeSetLinearSolver()")) {return 1;}

    /* Set the Jacobian routine */
    retval = CVodeSetJacFn(cvode_mem, linearosc_jac);
    if (check_int_status(retval, "CVodeSetJacFn()")) {return 1;}

    sunrealtype t1 = SUN_RCONST(3.1415926535897932384626433);

    retval = CVodeSetStopTime(cvode_mem, t1);
    if (check_int_status(retval, "CVodeSetStopTime()")) {return 1;}

    while (t < t1) {
        /* Advance the solution. */
        retval = CVode(cvode_mem, t1, y, &t, CV_ONE_STEP);
        if (retval != CV_SUCCESS && retval != CV_TSTOP_RETURN) {
            fprintf(stderr, "retval=%d\n", retval);
            retval = -1;
            break;
        }
        else {
            retval = 0;
        }
    }

    double xfinal = NV_Ith_S(y, 0);
    double yfinal = NV_Ith_S(y, 1);

    /* Free memory */
    N_VDestroy(y);
    CVodeFree(&cvode_mem);
    SUNLinSolFree(LS);
    SUNMatDestroy(A);
    SUNContext_Free(&sunctx);

    if (retval) {
        /* Something went wrong while solving. */
        return retval;
    }
    /* Check that we got the expected final value of (-1, 0). */
    if (fabs(xfinal + 1.0) > 1e-12 || fabs(yfinal) > 1e-12) {
        return 1;
    }
    return 0;
}
