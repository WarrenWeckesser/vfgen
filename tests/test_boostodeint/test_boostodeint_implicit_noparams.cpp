
#include <cmath>
#include <boost/numeric/odeint.hpp>

#include "linearosc_vf.h"

using namespace std;
using namespace boost::numeric::odeint;


int main(int argc, char *argv[])
{
    state_type z(2);
    double t0 = 0.0;
    double tfinal = M_PI;
    double dt = 0.1;

    z[0] = 1.0;
    z[1] = 0.0;

    auto linearosc = linearosc_vf();
    auto sys = make_pair(linearosc,
                         [&linearosc](const state_type &x_, matrix_type &J_, const double &t_, state_type &dfdt_)
                             {linearosc.jac(x_, J_, t_, dfdt_);});

    size_t nsteps = integrate_adaptive(
        make_dense_output<rosenbrock4<double>>(1e-12, 1e-10),
        sys, z, t0, tfinal, dt);

    if ((fabs(z[0] + 1) < 1e-9) && (fabs(z[1]) < 1e-9)) {
        return 0;
    }
    else {
        return -1;
    }
}
