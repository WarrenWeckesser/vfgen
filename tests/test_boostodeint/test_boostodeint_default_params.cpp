
#include <cmath>
#include <boost/numeric/odeint.hpp>
#include <boost/numeric/odeint/stepper/bulirsch_stoer.hpp>

#include "linearoscp_vf.h"

using namespace std;
using namespace boost::numeric::odeint;


int main(int argc, char *argv[])
{
    state_type z{1.0, 0.0};
    double t0 = 0.0;
    double tfinal = M_PI;
    double dt0 = 0.1;

    auto linearoscp = linearoscp_vf(2.0);

    bulirsch_stoer<state_type> stepper(1e-12, 1e-10, 1.0, 1.0);
    size_t nsteps = integrate_adaptive(stepper, linearoscp, z, t0, tfinal, dt0);

    if ((fabs(z[0] - 1) < 1e-9) && (fabs(z[1]) < 1e-9)) {
        return 0;
    }
    else {
        return -1;
    }
}
