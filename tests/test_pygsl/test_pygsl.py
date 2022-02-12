#
# test_pygsl.py
#

import numpy as np
from numpy.testing import assert_allclose
from pygsl import odeiv
import linearosc


# Create the GSL ODE solver
ndim = 2
step    = odeiv.step_rk8pd(ndim, linearosc.vectorfield, linearosc.jacobian)
control = odeiv.control_y_new(step, eps_abs=1e-12, eps_rel=1e-9)
evolve  = odeiv.evolve(step, control, ndim)

stoptime = np.pi
# Initial step size is stoptime/500
h = stoptime/500.0
t = 0
y = [1.0, 0.0]
# Call evolve.apply(...) until the solution reaches stoptime
while t < stoptime:
    t, h, y = evolve.apply(t, stoptime, h, y)

assert_allclose(y, [-1, 0], rtol=1e-8, atol=1e-9)
