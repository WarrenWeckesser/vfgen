/*
 *  test_taylor.c
 */

#include <stdio.h>
#include "one_over_x_taylor4.h"

int main (int argc, char *argv[])
{
    double x[1] = {2.0};
    double derivs[4][1];
    double ref[4] = {0.5, -1.0/8, 3.0/32, -15.0/128};

    one_over_x_derivs4(derivs, x, NULL);

    for (int k = 0; k < 4; ++k) {
        if (derivs[k][0] != ref[k]) {
            fprintf(stderr, "Computed deriv does not equal reference value:\n");
            fprintf(stderr, "k = %d   derivs[%d][0] = %20.15f  ref[%d] = %20.15f\n",
                    k, k, derivs[k][0], k, ref[k]);
            return -1;
        }
    }
    return 0;
}
